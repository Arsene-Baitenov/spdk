#include "bdev_raid.h"

#include "spdk/env.h"
#include "spdk/thread.h"
#include "spdk/string.h"
#include "spdk/util.h"
#include "spdk/likely.h"
#include "spdk/log.h"
#include "spdk/accel.h"

#define RAID5_MAX_STRIPES 32
// #define RAID5_MAX_STRIPES 4

struct chunk {
	/* Corresponds to base_bdev index */
	uint8_t index;

	/* Array of iovecs */
	struct iovec *iovs;

	/* Number of used iovecs */
	int iovcnt;

	/* Total number of available iovecs in the array */
	int iovcnt_max;

	bool enabled;
	uint64_t offset_blocks;
	uint64_t num_blocks;
	struct iovec *req_iovs;
	int req_iovcnt;
};

struct stripe_request;
typedef void (*stripe_req_xor_cb)(struct stripe_request *stripe_req, int status);
typedef void (*stripe_req_cb)(struct stripe_request *stripe_req, enum spdk_bdev_io_status status);

enum stripe_request_type {
	READ,
	WRITE
};

struct stripe_request {
	enum stripe_request_type type;

	struct raid5_io_channel *r5ch;

	/* The associated raid_bdev_io */
	struct raid_bdev_io *raid_io;

	/* The stripe's index in the raid array. */
	uint64_t stripe_index;

	stripe_req_cb cb;

	/* extra iovec for calculating parity in read-modify-write case */
	struct iovec extra;

    struct {
        /* Array of buffers for reading chunk data */
        void **chunk_buffers;

        /* zero buffer to pad to full strip */
        void *zero_buffer;
    } buff;

    struct {
        /* Array of iovec iterators for each chunk */
        struct spdk_ioviter *chunk_iov_iters;

        /* For iterating over chunk iovecs during xor calculation */
        void **chunk_buffers;
        struct iovec **chunk_iovs;
        size_t *chunk_iovcnt;
		uint32_t count;

        size_t len;
		size_t remaining;
		int status;
		stripe_req_xor_cb cb;
    } xor;

	TAILQ_ENTRY(stripe_request) link;

	/* Array of chunks corresponding to base_bdevs */
	struct chunk chunks[0];
};

struct stripe_lock {
	struct spdk_spinlock lock;
	bool locked;
};

struct raid5_info {
	/* The parent raid bdev */
	struct raid_bdev *raid_bdev;

	/* Number of data blocks in a stripe (without parity) */
	uint64_t stripe_blocks;

	/* Number of stripes on this array */
	uint64_t total_stripes;

	/* Alignment for buffer allocation */
	size_t buf_alignment;

	/* Array of locks for every stripe */
	struct stripe_lock *stripe_locks;
};

struct raid5_io_channel {
	/* All available stripe requests on this channel */
    TAILQ_HEAD(, stripe_request) free_stripe_requests;

	/* free_stripe_requests queue lock */
	struct spdk_spinlock lock;

	/* accel_fw channel */
	struct spdk_io_channel *accel_ch;
};

// various auxiliary functions

#define __CHUNK_IN_RANGE(req, c) \
	c < req->chunks + raid5_ch_to_r5_info(req->r5ch)->raid_bdev->num_base_bdevs

#define FOR_EACH_CHUNK_FROM(req, c, from) \
	for (c = from; __CHUNK_IN_RANGE(req, c); c++)

#define FOR_EACH_CHUNK(req, c) \
	FOR_EACH_CHUNK_FROM(req, c, req->chunks)

static inline struct raid5_info *
raid5_ch_to_r5_info(struct raid5_io_channel *r5ch)
{
	return spdk_io_channel_get_io_device(spdk_io_channel_from_ctx(r5ch));
}

static inline struct stripe_request *
raid5_chunk_stripe_req(struct chunk *chunk)
{
	return SPDK_CONTAINEROF((chunk - chunk->index), struct stripe_request, chunks);
}

static inline uint8_t
raid5_stripe_data_chunks_num(struct raid_bdev *raid_bdev)
{
	return raid_bdev->min_base_bdevs_operational;
}

static inline uint8_t
raid5_parity_chunk_idx(struct raid_bdev *raid_bdev, uint64_t stripe_index)
{
	return raid5_stripe_data_chunks_num(raid_bdev) - stripe_index % raid_bdev->num_base_bdevs;
}

static struct stripe_request *
raid5_get_stripe_request(struct raid5_io_channel *r5ch)
{
	struct stripe_request *stripe_req = NULL;

	spdk_spin_lock(&r5ch->lock);
	stripe_req = TAILQ_FIRST(&r5ch->free_stripe_requests);
	if (stripe_req) {
		TAILQ_REMOVE(&r5ch->free_stripe_requests, stripe_req, link);
	}
	spdk_spin_unlock(&r5ch->lock);

	return stripe_req;
}

static inline void
raid5_stripe_request_release(struct stripe_request *stripe_req)
{
	struct raid5_io_channel *r5ch = stripe_req->r5ch;

	spdk_spin_lock(&r5ch->lock);
    TAILQ_INSERT_HEAD(&r5ch->free_stripe_requests, stripe_req, link);
	spdk_spin_unlock(&r5ch->lock);
}

static inline uint8_t
raid5_start_chunk_idx(struct stripe_request *stripe_req)
{
	uint8_t start_idx = (stripe_req->raid_io->offset_blocks >>
		stripe_req->raid_io->raid_bdev->strip_size_shift) %
		raid5_stripe_data_chunks_num(stripe_req->raid_io->raid_bdev);
	if (raid5_parity_chunk_idx(stripe_req->raid_io->raid_bdev, stripe_req->stripe_index) <= start_idx) {
		++start_idx;
	}
	return start_idx;
}

static inline uint8_t
raid5_end_chunk_idx(struct stripe_request *stripe_req)
{
	uint8_t end_idx = ((stripe_req->raid_io->offset_blocks + stripe_req->raid_io->num_blocks - 1) >>
		stripe_req->raid_io->raid_bdev->strip_size_shift) %
		raid5_stripe_data_chunks_num(stripe_req->raid_io->raid_bdev);
	if (raid5_parity_chunk_idx(stripe_req->raid_io->raid_bdev, stripe_req->stripe_index) <= end_idx) {
		++end_idx;
	}
	return end_idx;
}

static inline uint64_t
raid5_chunk_offb(struct stripe_request *stripe_req, uint8_t idx)
{
	return idx == raid5_start_chunk_idx(stripe_req) ?
			(stripe_req->raid_io->offset_blocks &
			(stripe_req->raid_io->raid_bdev->strip_size - 1)) : 0;
}

static inline uint64_t
raid5_chunk_endb(struct stripe_request *stripe_req, uint8_t idx)
{
	return idx == raid5_end_chunk_idx(stripe_req) ?
			(stripe_req->raid_io->num_blocks +
			(stripe_req->raid_io->offset_blocks &
			(stripe_req->raid_io->raid_bdev->strip_size - 1)) +
			stripe_req->raid_io->raid_bdev->strip_size - 1)
				% stripe_req->raid_io->raid_bdev->strip_size + 1 :
				stripe_req->raid_io->raid_bdev->strip_size;
}

static void
raid5_chunks_disable(struct stripe_request *stripe_req)
{
	struct chunk *chunk;
	FOR_EACH_CHUNK(stripe_req, chunk) {
		chunk->enabled = false;
	}
}

static enum spdk_bdev_io_status
raid5_return_code_to_io_status(int ret)
{
	if (ret == 0) {
		return SPDK_BDEV_IO_STATUS_SUCCESS;
	}
	if (ret == -ENOMEM) {
		return SPDK_BDEV_IO_STATUS_NOMEM;
	}
	return SPDK_BDEV_IO_STATUS_FAILED;
}

// lock/unlock stripe functions

static inline int
raid5_try_lock_stripe(struct raid5_info *r5_info, uint64_t stripe_idx)
{
	int ret = 0;

	spdk_spin_lock(&r5_info->stripe_locks[stripe_idx].lock);
	if (r5_info->stripe_locks[stripe_idx].locked) {
		ret = -ENOMEM;
	} else {
		r5_info->stripe_locks[stripe_idx].locked = true;
	}
	spdk_spin_unlock(&r5_info->stripe_locks[stripe_idx].lock);
	
	return ret;
}

static inline void
raid5_unlock_stripe(struct raid5_info *r5_info, uint64_t stripe_idx)
{
	r5_info->stripe_locks[stripe_idx].locked = false;
}

// xor processing functions

static void
raid5_xor_done(struct stripe_request *stripe_req)
{
	if (stripe_req->xor.status != 0) {
		SPDK_ERRLOG("Xor failed with error: %s\n", spdk_strerror(-stripe_req->xor.status));
	}

	stripe_req->xor.cb(stripe_req, stripe_req->xor.status);
}

static void raid5_xor_continue(struct stripe_request *stripe_req);

static void
raid5_xor_done_part(void *_stripe_req, int status)
{
	struct stripe_request *stripe_req = _stripe_req;

	assert(stripe_req->xor.remaining >= stripe_req->xor.len);
	stripe_req->xor.remaining -= stripe_req->xor.len;

	if (status != 0) {
		stripe_req->xor.status = status;
	}

	if (stripe_req->xor.remaining > 0 && stripe_req->xor.status == 0) {
		stripe_req->xor.len = spdk_ioviter_nextv(
			stripe_req->xor.chunk_iov_iters,
			stripe_req->xor.chunk_buffers);
		if (stripe_req->xor.len == 0) {
			stripe_req->xor.status = -EINVAL;
			raid5_xor_done(stripe_req);
			return;
		}
		raid5_xor_continue(stripe_req);
	} else {
		raid5_xor_done(stripe_req);
	}
}

static void
raid5_xor_continue(struct stripe_request *stripe_req)
{
	struct raid5_io_channel *r5ch = stripe_req->r5ch;
	uint32_t nsrsc = stripe_req->xor.count - 1;
	int ret;

	assert(stripe_req->xor.len > 0);

	ret = spdk_accel_submit_xor(
		r5ch->accel_ch,
		stripe_req->xor.chunk_buffers[nsrsc],
		stripe_req->xor.chunk_buffers,
		nsrsc,
		stripe_req->xor.len,
		raid5_xor_done_part,
		stripe_req);

	if (spdk_unlikely(ret)) {
		stripe_req->xor.status = ret;
		raid5_xor_done(stripe_req);
	}
}

static void
raid5_xor_start(struct stripe_request *stripe_req, stripe_req_xor_cb cb)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;

	stripe_req->xor.len = spdk_ioviter_firstv(
		stripe_req->xor.chunk_iov_iters,
		stripe_req->xor.count,
		stripe_req->xor.chunk_iovs,
		stripe_req->xor.chunk_iovcnt,
		stripe_req->xor.chunk_buffers);

	stripe_req->xor.remaining = raid_bdev->strip_size * raid_bdev->bdev.blocklen;
	stripe_req->xor.status = 0;
	stripe_req->xor.cb = cb;

	raid5_xor_continue(stripe_req);
}

// stripe request completion functions

static void
raid5_stripe_request_complete(
	struct stripe_request *stripe_req,
	enum spdk_bdev_io_status status)
{
	if (stripe_req->cb != NULL) {
		stripe_req->cb(stripe_req, status);
	} else {
		struct raid_bdev_io *raid_io = stripe_req->raid_io;

		raid5_stripe_request_release(stripe_req);

		raid_bdev_io_complete(raid_io, status);
	}
}

static bool
raid5_stripe_request_complete_part(
	struct stripe_request *stripe_req,
	uint64_t completed,
	enum spdk_bdev_io_status status)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;

	assert(raid_io->base_bdev_io_remaining >= completed);
	raid_io->base_bdev_io_remaining -= completed;

	if (status != SPDK_BDEV_IO_STATUS_SUCCESS) {
		raid_io->base_bdev_io_status = status;
	}

	if (raid_io->base_bdev_io_remaining == 0) {
		raid5_stripe_request_complete(stripe_req, raid_io->base_bdev_io_status);
		return true;
	} else {
		return false;
	}
}

// raid5 iovecs mapping helper functions

static int
raid5_chunk_set_iovcnt(struct chunk *chunk, int iovcnt)
{
	if (iovcnt > chunk->iovcnt_max) {
		struct iovec *iovs = chunk->iovs;

		iovs = realloc(iovs, iovcnt * sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}
		chunk->iovs = iovs;
		chunk->iovcnt_max = iovcnt;
	}
	chunk->iovcnt = iovcnt;

	return 0;
}

static int
raid5_map_iovecs(struct stripe_request *stripe_req, struct iovec *iovs, int iovcnt)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct chunk *chunk;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	int iov_idx = 0;
	uint64_t offset = 0;
	uint64_t iov_offset = 0;

	for (int i=start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		chunk = &stripe_req->chunks[i];
		int chunk_iovcnt = 0;
		uint64_t len = raid_bdev->strip_size * raid_bdev->bdev.blocklen;
		uint64_t off = iov_offset;
		int ret = 0;

		for (int i = iov_idx; i < iovcnt; i++) {
			chunk_iovcnt++;
			off += iovs[i].iov_len;
			if (off >= offset + len) {
				break;
			}
		}

		ret = raid5_chunk_set_iovcnt(chunk, chunk_iovcnt);
		if (ret) {
			return ret;
		}

		for (int i = 0; i < chunk_iovcnt; i++) {
			struct iovec *chunk_iov = &chunk->iovs[i];
			const struct iovec *iov = &iovs[iov_idx];
			size_t chunk_iov_offset = offset - iov_offset;

			chunk_iov->iov_base = iov->iov_base + chunk_iov_offset;
			chunk_iov->iov_len = spdk_min(len, iov->iov_len - chunk_iov_offset);
			offset += chunk_iov->iov_len;
			len -= chunk_iov->iov_len;

			if (offset >= iov_offset + iov->iov_len) {
				iov_idx++;
				iov_offset += iov->iov_len;
			}
		}

		if (spdk_unlikely(len > 0)) {
			return -EINVAL;
		}

		chunk->offset_blocks = 0;
		chunk->num_blocks = raid_bdev->strip_size;
		chunk->req_iovs = chunk->iovs;
		chunk->req_iovcnt = chunk->iovcnt;
		chunk->enabled = true;
	}
	return 0;
}

// functions of processing requests to base bdevs

static void
raid5_chunk_complete(struct chunk *chunk, enum spdk_bdev_io_status status)
{
	struct stripe_request *stripe_req = raid5_chunk_stripe_req(chunk);

	raid5_stripe_request_complete_part(stripe_req, 1, status);
}

static void
raid5_chunk_complete_bdev_io(struct spdk_bdev_io *bdev_io, bool success, void *cb_arg)
{
	struct chunk *chunk = cb_arg;
	enum spdk_bdev_io_status status = success ? SPDK_BDEV_IO_STATUS_SUCCESS :
					  SPDK_BDEV_IO_STATUS_FAILED;

	spdk_bdev_free_io(bdev_io);

	raid5_chunk_complete(chunk, status);
}

static void raid5_submit_chunks(struct stripe_request *stripe_req);

static void
raid5_submit_chunk_retry(void *_raid_io)
{
	struct raid_bdev_io *raid_io = _raid_io;
	struct stripe_request *stripe_req = raid_io->module_private;

	raid5_submit_chunks(stripe_req);
}

static inline void
raid5_init_ext_io_opts(struct spdk_bdev_ext_io_opts *opts, struct raid_bdev_io *raid_io)
{
	memset(opts, 0, sizeof(*opts));
	opts->size = sizeof(*opts);
	opts->memory_domain = raid_io->memory_domain;
	opts->memory_domain_ctx = raid_io->memory_domain_ctx;
	opts->metadata = raid_io->md_buf;
}

static int
raid5_submit_chunk(struct chunk *chunk)
{
	struct stripe_request *stripe_req = raid5_chunk_stripe_req(chunk);
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct raid_base_bdev_info *base_info = &raid_bdev->base_bdev_info[chunk->index];
	struct spdk_io_channel *base_ch = raid_bdev_channel_get_base_channel(raid_io->raid_ch, chunk->index);
	struct spdk_bdev_ext_io_opts io_opts;
	uint64_t offset_blocks;
	uint64_t num_blocks;
	int ret;

	if (!chunk->enabled) {
		raid5_chunk_complete(chunk, SPDK_BDEV_IO_STATUS_SUCCESS);
		return 0;
	}

	num_blocks = chunk->num_blocks;
	offset_blocks = (stripe_req->stripe_index << raid_bdev->strip_size_shift) + chunk->offset_blocks;

	raid5_init_ext_io_opts(&io_opts, raid_io);
	io_opts.metadata = NULL;
	
	switch (stripe_req->type) {
	case WRITE:
		ret = raid_bdev_writev_blocks_ext(
			base_info,
			base_ch,
			chunk->req_iovs,
			chunk->req_iovcnt,
			offset_blocks,
			num_blocks,
			raid5_chunk_complete_bdev_io,
			chunk,
			&io_opts);
		break;
	case READ:
		ret = raid_bdev_readv_blocks_ext(
			base_info,
			base_ch,
			chunk->req_iovs,
			chunk->req_iovcnt,
			offset_blocks,
			num_blocks,
			raid5_chunk_complete_bdev_io,
			chunk,
			&io_opts);
		break;
	default:
		assert(false);
		ret = -EINVAL;
		break;
	}

	if (spdk_unlikely(ret)) {
		if (ret == -ENOMEM) {
			raid_bdev_queue_io_wait(
				raid_io,
				spdk_bdev_desc_get_bdev(base_info->desc),
				base_ch,
				raid5_submit_chunk_retry);
		} else {
			uint64_t not_submitted = raid_bdev->num_base_bdevs -
									raid_io->base_bdev_io_submitted;

			raid5_stripe_request_complete_part(
				stripe_req,
				not_submitted,
				SPDK_BDEV_IO_STATUS_FAILED);
		}
	}

	return ret;
}

static void
raid5_submit_chunks(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct chunk *start = &stripe_req->chunks[raid_io->base_bdev_io_submitted];
	struct chunk *chunk;

	FOR_EACH_CHUNK_FROM(stripe_req, chunk, start) {
		if (spdk_likely(raid5_submit_chunk(chunk) == 0)) {
			raid_io->base_bdev_io_submitted++;
		} else {
			break;
		}
	}
}

static void
raid5_submit_stripe_request(
	struct stripe_request *stripe_req,
	enum stripe_request_type type,
	stripe_req_cb cb)
{
	stripe_req->type = type;
	stripe_req->raid_io->module_private = stripe_req;
	stripe_req->raid_io->base_bdev_io_remaining = stripe_req->raid_io->raid_bdev->num_base_bdevs;
	stripe_req->raid_io->base_bdev_io_submitted = 0;
	stripe_req->cb = cb;

	raid5_submit_chunks(stripe_req);
}

// write request processing function

static int
raid5_remap_iovecs_rmw(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct iovec *iovs;
	int iovcnt;
	int start_iov_idx_ofs = 0;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint64_t strip_size_bytes = raid_bdev->strip_size * raid_bdev->bdev.blocklen;
	uint8_t c = 0;
	int ret = 0;

	raid5_chunks_disable(stripe_req);

	{
		iovcnt = raid_io->iovcnt + (start_offb != 0 ? 1 : 0) + (end_endb != raid_bdev->strip_size ? 1 : 0);

		iovs = calloc(iovcnt, sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}

		if (start_offb != 0) {
			start_iov_idx_ofs = 1;
			iovs[0].iov_base = stripe_req->buff.zero_buffer;
			iovs[0].iov_len = start_offb * raid_bdev->bdev.blocklen;
		}
		if (end_endb != raid_bdev->strip_size) {
			iovs[iovcnt - 1].iov_base = stripe_req->buff.zero_buffer + (end_endb * raid_bdev->bdev.blocklen);
			iovs[iovcnt - 1].iov_len = (raid_bdev->strip_size - end_endb) * raid_bdev->bdev.blocklen;
		}

		for (int i = 0; i < raid_io->iovcnt; ++i) {
			iovs[i + start_iov_idx_ofs].iov_base = raid_io->iovs[i].iov_base;
			iovs[i + start_iov_idx_ofs].iov_len = raid_io->iovs[i].iov_len;
		}
	}

	ret = raid5_map_iovecs(stripe_req, iovs, iovcnt);
	if (ret) {
		free(iovs);
		return ret;
	}
	free(iovs);

	if (start_offb != 0) {
		stripe_req->chunks[start].offset_blocks = start_offb;
		stripe_req->chunks[start].num_blocks -= start_offb;
		stripe_req->chunks[start].req_iovs = &stripe_req->chunks[start].iovs[1];
		stripe_req->chunks[start].req_iovcnt--;
	}
	if (end_endb != raid_bdev->strip_size) {
		stripe_req->chunks[end].num_blocks -= end_endb;
		stripe_req->chunks[end].req_iovcnt--;
	}

	stripe_req->chunks[parity].iovs[0].iov_base = stripe_req->buff.chunk_buffers[parity];
	stripe_req->chunks[parity].iovs[0].iov_len = strip_size_bytes;
	stripe_req->chunks[parity].iovcnt = 1;
	stripe_req->chunks[parity].offset_blocks = 0;
	stripe_req->chunks[parity].num_blocks = raid_bdev->strip_size;
	stripe_req->chunks[parity].req_iovs = stripe_req->chunks[parity].iovs;
	stripe_req->chunks[parity].req_iovcnt = stripe_req->chunks[parity].iovcnt;
	stripe_req->chunks[parity].enabled = true;

	stripe_req->xor.chunk_iovs[c] = &stripe_req->extra;
	stripe_req->xor.chunk_iovcnt[c] = 1;
	++c;

	for (int i = start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[i].iovs;
		stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[i].iovcnt;
		++c;
	}
	stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[parity].iovs;
	stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[parity].iovcnt;

	stripe_req->xor.count = c + 1;

	return 0;
}

static int
raid5_prepare_xor_rmw(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint8_t c = 0;

	raid5_chunks_disable(stripe_req);

	for (int i = start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[i].iovs;
		stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[i].iovcnt;
		++c;
	}
	stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[parity].iovs;
	stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[parity].iovcnt;
	++c;

	stripe_req->xor.chunk_iovs[c] = &stripe_req->extra;
	stripe_req->xor.chunk_iovcnt[c] = 1;

	stripe_req->xor.count = c + 1;

	return 0;
}

static int
raid5_map_iovecs_rmw(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint64_t strip_size_bytes = raid_bdev->strip_size * raid_bdev->bdev.blocklen;

	for (int i = start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		stripe_req->chunks[i].iovs[0].iov_base = stripe_req->buff.chunk_buffers[i];
		stripe_req->chunks[i].iovs[0].iov_len = strip_size_bytes;
		stripe_req->chunks[i].iovcnt = 1;
		stripe_req->chunks[i].offset_blocks = 0;
		stripe_req->chunks[i].num_blocks = raid_bdev->strip_size;
		stripe_req->chunks[i].req_iovs = stripe_req->chunks[i].iovs;
		stripe_req->chunks[i].req_iovcnt = stripe_req->chunks[i].iovcnt;
		stripe_req->chunks[i].enabled = true;
	}

	stripe_req->chunks[parity].iovs[0].iov_base = stripe_req->buff.chunk_buffers[parity];
	stripe_req->chunks[parity].iovs[0].iov_len = strip_size_bytes;
	stripe_req->chunks[parity].iovcnt = 1;
	stripe_req->chunks[parity].offset_blocks = 0;
	stripe_req->chunks[parity].num_blocks = raid_bdev->strip_size;
	stripe_req->chunks[parity].req_iovs = stripe_req->chunks[parity].iovs;
	stripe_req->chunks[parity].req_iovcnt = stripe_req->chunks[parity].iovcnt;
	stripe_req->chunks[parity].enabled = true;
	return 0;
}

static int
raid5_remap_iovecs_write_reconstruct(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct iovec *iovs;
	int iovcnt;
	int start_iov_idx_ofs = 0;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint64_t strip_size_bytes = raid_bdev->strip_size * raid_bdev->bdev.blocklen;
	uint8_t c = 0;
	int ret = 0;

	raid5_chunks_disable(stripe_req);

	{
		iovcnt = raid_io->iovcnt + (start_offb != 0 ? 1 : 0) + (end_endb != raid_bdev->strip_size ? 1 : 0);

		iovs = calloc(iovcnt, sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}

		if (start_offb != 0) {
			start_iov_idx_ofs = 1;
			iovs[0].iov_base = stripe_req->buff.chunk_buffers[start];
			iovs[0].iov_len = start_offb * raid_bdev->bdev.blocklen;
		}
		if (end_endb != raid_bdev->strip_size) {
			iovs[iovcnt - 1].iov_base = stripe_req->buff.chunk_buffers[end] + (end_endb * raid_bdev->bdev.blocklen);
			iovs[iovcnt - 1].iov_len = (raid_bdev->strip_size - end_endb) * raid_bdev->bdev.blocklen;
		}

		for (int i = 0; i < raid_io->iovcnt; ++i) {
			iovs[i + start_iov_idx_ofs].iov_base = raid_io->iovs[i].iov_base;
			iovs[i + start_iov_idx_ofs].iov_len = raid_io->iovs[i].iov_len;
		}
	}

	ret = raid5_map_iovecs(stripe_req, iovs, iovcnt);
	if (ret) {
		free(iovs);
		return ret;
	}
	free(iovs);

	for (int i = start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		if (raid_bdev_channel_get_base_channel(raid_io->raid_ch, i) == NULL) {
			stripe_req->chunks[i].enabled = false;
			break;
		}
	}

	stripe_req->chunks[parity].iovs[0].iov_base = stripe_req->buff.chunk_buffers[parity];
	stripe_req->chunks[parity].iovs[0].iov_len = strip_size_bytes;
	stripe_req->chunks[parity].iovcnt = 1;
	stripe_req->chunks[parity].offset_blocks = 0;
	stripe_req->chunks[parity].num_blocks = raid_bdev->strip_size;
	stripe_req->chunks[parity].req_iovs = stripe_req->chunks[parity].iovs;
	stripe_req->chunks[parity].req_iovcnt = stripe_req->chunks[parity].iovcnt;
	stripe_req->chunks[parity].enabled = true;

	for (int i = 0; i < raid_bdev->num_base_bdevs; i++) {
		if (i == parity) {
			continue;
		}
		stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[i].iovs;
		stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[i].iovcnt;
		++c;
	}
	stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[parity].iovs;
	stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[parity].iovcnt;
	
	stripe_req->xor.count = c + 1;
	return 0;
}

static int
raid5_map_iovecs_write_reconstruct(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint64_t strip_size_bytes = raid_bdev->strip_size * raid_bdev->bdev.blocklen;

	{
		if (start_offb != 0) {
			stripe_req->chunks[start].iovs[0].iov_base = stripe_req->buff.chunk_buffers[start];
			stripe_req->chunks[start].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[start].iovcnt = 1;
			stripe_req->chunks[start].offset_blocks = 0;
			stripe_req->chunks[start].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[start].req_iovs = stripe_req->chunks[start].iovs;
			stripe_req->chunks[start].req_iovcnt = stripe_req->chunks[start].iovcnt;
			stripe_req->chunks[start].enabled = true;
		}
		if (end_endb != raid_bdev->strip_size) {
			stripe_req->chunks[end].iovs[0].iov_base = stripe_req->buff.chunk_buffers[end];
			stripe_req->chunks[end].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[end].iovcnt = 1;
			stripe_req->chunks[end].offset_blocks = 0;
			stripe_req->chunks[end].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[end].req_iovs = stripe_req->chunks[end].iovs;
			stripe_req->chunks[end].req_iovcnt = stripe_req->chunks[end].iovcnt;
			stripe_req->chunks[end].enabled = true;
		}
	}

	{
		for (int i = 0; i < start; i = (i + 1 == parity ? i + 2 : i + 1)) {
			stripe_req->chunks[i].iovs[0].iov_base = stripe_req->buff.chunk_buffers[i];
			stripe_req->chunks[i].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[i].iovcnt = 1;
			stripe_req->chunks[i].offset_blocks = 0;
			stripe_req->chunks[i].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[i].req_iovs = stripe_req->chunks[i].iovs;
			stripe_req->chunks[i].req_iovcnt = stripe_req->chunks[i].iovcnt;
			stripe_req->chunks[i].enabled = true;
		}

		for (int i = end + 1; i < raid_bdev->num_base_bdevs; i = (i + 1 == parity ? i + 2 : i + 1)) {
			stripe_req->chunks[i].iovs[0].iov_base = stripe_req->buff.chunk_buffers[i];
			stripe_req->chunks[i].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[i].iovcnt = 1;
			stripe_req->chunks[i].offset_blocks = 0;
			stripe_req->chunks[i].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[i].req_iovs = stripe_req->chunks[i].iovs;
			stripe_req->chunks[i].req_iovcnt = stripe_req->chunks[i].iovcnt;
			stripe_req->chunks[i].enabled = true;
		}
	}
	return 0;
}

static int
raid5_map_iovecs_write(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct iovec *iovs;
	int iovcnt;
	int start_iov_idx_ofs = 0;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	int ret = 0;

	{
		iovcnt = raid_io->iovcnt + (start_offb != 0 ? 1 : 0) + (end_endb != raid_bdev->strip_size ? 1 : 0);

		iovs = calloc(iovcnt, sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}

		if (start_offb != 0) {
			start_iov_idx_ofs = 1;
			iovs[0].iov_base = stripe_req->buff.chunk_buffers[start];
			iovs[0].iov_len = start_offb * raid_bdev->bdev.blocklen;
		}
		if (end_endb != raid_bdev->strip_size) {
			iovs[iovcnt - 1].iov_base = stripe_req->buff.chunk_buffers[end] + (end_endb * raid_bdev->bdev.blocklen);
			iovs[iovcnt - 1].iov_len = (raid_bdev->strip_size - end_endb) * raid_bdev->bdev.blocklen;
		}

		for (int i = 0; i < raid_io->iovcnt; ++i) {
			iovs[i + start_iov_idx_ofs].iov_base = raid_io->iovs[i].iov_base;
			iovs[i + start_iov_idx_ofs].iov_len = raid_io->iovs[i].iov_len;
		}
	}

	ret = raid5_map_iovecs(stripe_req, iovs, iovcnt);
	if (ret) {
		free(iovs);
		return ret;
	}

	if (start_offb != 0) {
		stripe_req->chunks[start].offset_blocks = start_offb;
		stripe_req->chunks[start].num_blocks -= start_offb;
		stripe_req->chunks[start].req_iovs = &stripe_req->chunks[start].iovs[1];
		stripe_req->chunks[start].req_iovcnt--;
	}
	if (end_endb != raid_bdev->strip_size) {
		stripe_req->chunks[end].num_blocks -= end_endb;
		stripe_req->chunks[end].req_iovcnt--;
	}

	free(iovs);
	return 0;
}

static void
raid5_write_unlock_and_complete(
	struct stripe_request *stripe_req,
	enum spdk_bdev_io_status status)
{
	struct raid5_info *r5_info = stripe_req->raid_io->raid_bdev->module_private;

	stripe_req->cb = NULL;
	raid5_unlock_stripe(r5_info, stripe_req->stripe_index);
	raid5_stripe_request_complete(stripe_req, status);
}

static void
raid5_rmw_write_xor_cb(
	struct stripe_request *stripe_req,
	int status)
{
	if (status != 0) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(status));
		return;
	}

	raid5_submit_stripe_request(stripe_req, WRITE, raid5_write_unlock_and_complete);
}

static void
raid5_rmw_read_xor_cb(
	struct stripe_request *stripe_req,
	int status)
{
	int ret = 0;

	if (status != 0) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(status));
		return;
	}

	ret = raid5_remap_iovecs_rmw(stripe_req);
	if (ret) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(ret));
		return;
	}

	raid5_xor_start(stripe_req, raid5_rmw_write_xor_cb);
}

static void
raid5_rmw_cb(
	struct stripe_request *stripe_req,
	enum spdk_bdev_io_status status)
{
	int ret = 0;

	stripe_req->cb = NULL;

	if (status != SPDK_BDEV_IO_STATUS_SUCCESS) {
		raid5_write_unlock_and_complete(stripe_req, status);
		return;
	}

	ret = raid5_prepare_xor_rmw(stripe_req);
	if (ret) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(ret));
		return;
	}

	raid5_xor_start(stripe_req, raid5_rmw_read_xor_cb);
}

static void
raid5_write_reconstruct_xor_cb(
	struct stripe_request *stripe_req,
	int status)
{
	if (status != 0) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(status));
		return;
	}

	raid5_submit_stripe_request(stripe_req, WRITE, raid5_write_unlock_and_complete);
}

static void
raid5_write_reconstruct_cb(
	struct stripe_request *stripe_req,
	enum spdk_bdev_io_status status)
{
	int ret = 0;

	stripe_req->cb = NULL;

	if (status != SPDK_BDEV_IO_STATUS_SUCCESS) {
		raid5_write_unlock_and_complete(stripe_req, status);
		return;
	}

	ret = raid5_remap_iovecs_write_reconstruct(stripe_req);
	if (ret) {
		raid5_write_unlock_and_complete(
			stripe_req,
			raid5_return_code_to_io_status(ret));
		return;
	}

	raid5_xor_start(stripe_req, raid5_write_reconstruct_xor_cb);
}

static int
raid5_submit_write_request(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct raid5_info *r5_info = raid_bdev->module_private;
	struct spdk_io_channel *base_ch;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint8_t broken = raid_bdev->num_base_bdevs;
	stripe_req_cb cb = NULL;
	enum stripe_request_type type;
	int ret = 0;

	for (int i=start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		base_ch = raid_bdev_channel_get_base_channel(raid_io->raid_ch, i);
		if (base_ch == NULL) {
			broken = i;
			break;
		}
	}

	if (broken == raid_bdev->num_base_bdevs) {
		base_ch = raid_bdev_channel_get_base_channel(raid_io->raid_ch, parity);
		if (base_ch == NULL) {
			// broken parity strip (write)
			ret = raid5_map_iovecs_write(stripe_req);
			if (ret) {
				return ret;
			}
			cb = NULL;
			type = WRITE;
		} else {
			// read-modify-write (read_modify_write)
			ret = raid5_try_lock_stripe(r5_info, stripe_req->stripe_index);
			if (ret) {
				return ret;
			}
			ret = raid5_map_iovecs_rmw(stripe_req);
			if (ret) {
				raid5_unlock_stripe(r5_info, stripe_req->stripe_index);
				return ret;
			}
			cb = raid5_rmw_cb;
			type = READ;
		}
	} else {
		// broken req strip (write_reconstruct)
		ret = raid5_try_lock_stripe(r5_info, stripe_req->stripe_index);
		if (ret) {
			return ret;
		}
		ret = raid5_map_iovecs_write_reconstruct(stripe_req);
		if (ret) {
			raid5_unlock_stripe(r5_info, stripe_req->stripe_index);
			return ret;
		}
		cb = raid5_write_reconstruct_cb;
		type = READ;
	}

	raid5_submit_stripe_request(stripe_req, type, cb);

	return 0;
}

// read request processing function

static int
raid5_prepare_xor_read_reconstruct(struct stripe_request *stripe_req)
{
	struct raid_bdev *raid_bdev = stripe_req->raid_io->raid_bdev;
	uint8_t broken = 0;
	uint8_t c = 0;

	for (int i = 0; i < raid_bdev->num_base_bdevs; i++) {
		if (!stripe_req->chunks[i].enabled) {
			broken = i;
			continue;
		}
		stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[i].iovs;
		stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[i].iovcnt;
		++c;
	}
	stripe_req->xor.chunk_iovs[c] = stripe_req->chunks[broken].iovs;
	stripe_req->xor.chunk_iovcnt[c] = stripe_req->chunks[broken].iovcnt;
	
	stripe_req->xor.count = c + 1;
	return 0;
}

static int
raid5_map_iovecs_read_reconstruct(
	struct stripe_request *stripe_req,
	uint8_t broken)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct iovec *iovs;
	int iovcnt;
	int start_iov_idx_ofs = 0;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint64_t strip_size_bytes = raid_bdev->strip_size * raid_bdev->bdev.blocklen;
	int ret = 0;

	{
		iovcnt = raid_io->iovcnt + (start_offb != 0 ? 1 : 0) + (end_endb != raid_bdev->strip_size ? 1 : 0);

		iovs = calloc(iovcnt, sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}

		if (start_offb != 0) {
			start_iov_idx_ofs = 1;
			iovs[0].iov_base = stripe_req->buff.chunk_buffers[start];
			iovs[0].iov_len = start_offb * raid_bdev->bdev.blocklen;
		}
		if (end_endb != raid_bdev->strip_size) {
			iovs[iovcnt - 1].iov_base = stripe_req->buff.chunk_buffers[end] + (end_endb * raid_bdev->bdev.blocklen);
			iovs[iovcnt - 1].iov_len = (raid_bdev->strip_size - end_endb) * raid_bdev->bdev.blocklen;
		}

		for (int i = 0; i < raid_io->iovcnt; ++i) {
			iovs[i + start_iov_idx_ofs].iov_base = raid_io->iovs[i].iov_base;
			iovs[i + start_iov_idx_ofs].iov_len = raid_io->iovs[i].iov_len;
		}
	}

	ret = raid5_map_iovecs(stripe_req, iovs, iovcnt);
	if (ret) {
		free(iovs);
		return ret;
	}
	free(iovs);

	{
		for (int i = 0; i < start; i = (i + 1 == parity ? i + 2 : i + 1)) {
			stripe_req->chunks[i].iovs[0].iov_base = stripe_req->buff.chunk_buffers[i];
			stripe_req->chunks[i].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[i].iovcnt = 1;
			stripe_req->chunks[i].offset_blocks = 0;
			stripe_req->chunks[i].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[i].req_iovs = stripe_req->chunks[i].iovs;
			stripe_req->chunks[i].req_iovcnt = stripe_req->chunks[i].iovcnt;
			stripe_req->chunks[i].enabled = true;
		}

		for (int i = end + 1; i < raid_bdev->num_base_bdevs; i = (i + 1 == parity ? i + 2 : i + 1)) {
			stripe_req->chunks[i].iovs[0].iov_base = stripe_req->buff.chunk_buffers[i];
			stripe_req->chunks[i].iovs[0].iov_len = strip_size_bytes;
			stripe_req->chunks[i].iovcnt = 1;
			stripe_req->chunks[i].offset_blocks = 0;
			stripe_req->chunks[i].num_blocks = raid_bdev->strip_size;
			stripe_req->chunks[i].req_iovs = stripe_req->chunks[i].iovs;
			stripe_req->chunks[i].req_iovcnt = stripe_req->chunks[i].iovcnt;
			stripe_req->chunks[i].enabled = true;
		}
	}

	stripe_req->chunks[parity].iovs[0].iov_base = stripe_req->buff.chunk_buffers[parity];
	stripe_req->chunks[parity].iovs[0].iov_len = strip_size_bytes;
	stripe_req->chunks[parity].iovcnt = 1;
	stripe_req->chunks[parity].offset_blocks = 0;
	stripe_req->chunks[parity].num_blocks = raid_bdev->strip_size;
	stripe_req->chunks[parity].req_iovs = stripe_req->chunks[parity].iovs;
	stripe_req->chunks[parity].req_iovcnt = stripe_req->chunks[parity].iovcnt;
	stripe_req->chunks[parity].enabled = true;

	stripe_req->chunks[broken].enabled = false;

	return 0;
}

static void
raid5_read_reconstruct_xor_cb(
	struct stripe_request *stripe_req,
	int status)
{
	raid5_stripe_request_complete(
		stripe_req,
		raid5_return_code_to_io_status(status));
}

static void
raid5_read_reconstruct_cb(
	struct stripe_request *stripe_req,
	enum spdk_bdev_io_status status)
{
	int ret = 0;

	stripe_req->cb = NULL;

	if (status != SPDK_BDEV_IO_STATUS_SUCCESS) {
		raid5_stripe_request_complete(stripe_req, status);
		return;
	}

	ret = raid5_prepare_xor_read_reconstruct(stripe_req);
	if (ret) {
		raid5_stripe_request_complete(
			stripe_req,
			raid5_return_code_to_io_status(ret));
		return;
	}

	raid5_xor_start(stripe_req, raid5_read_reconstruct_xor_cb);
}

static int
raid5_map_iovecs_read(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct iovec *iovs;
	int iovcnt;
	int start_iov_idx_ofs = 0;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint64_t start_offb = raid5_chunk_offb(stripe_req, start);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint64_t end_endb = raid5_chunk_endb(stripe_req, end);
	int ret = 0;

	{
		iovcnt = raid_io->iovcnt + (start_offb != 0 ? 1 : 0) + (end_endb != raid_bdev->strip_size ? 1 : 0);

		iovs = calloc(iovcnt, sizeof(*iovs));
		if (!iovs) {
			return -ENOMEM;
		}

		if (start_offb != 0) {
			start_iov_idx_ofs = 1;
			iovs[0].iov_base = stripe_req->buff.chunk_buffers[start];
			iovs[0].iov_len = start_offb * raid_bdev->bdev.blocklen;
		}
		if (end_endb != raid_bdev->strip_size) {
			iovs[iovcnt - 1].iov_base = stripe_req->buff.chunk_buffers[end] + (end_endb * raid_bdev->bdev.blocklen);
			iovs[iovcnt - 1].iov_len = (raid_bdev->strip_size - end_endb) * raid_bdev->bdev.blocklen;
		}

		for (int i = 0; i < raid_io->iovcnt; ++i) {
			iovs[i + start_iov_idx_ofs].iov_base = raid_io->iovs[i].iov_base;
			iovs[i + start_iov_idx_ofs].iov_len = raid_io->iovs[i].iov_len;
		}
	}

	ret = raid5_map_iovecs(stripe_req, iovs, iovcnt);
	if (ret) {
		free(iovs);
		return ret;
	}

	free(iovs);
	return 0;
}

static int
raid5_submit_read_request(struct stripe_request *stripe_req)
{
	struct raid_bdev_io *raid_io = stripe_req->raid_io;
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct spdk_io_channel *base_ch;
	uint8_t start = raid5_start_chunk_idx(stripe_req);
	uint8_t end = raid5_end_chunk_idx(stripe_req);
	uint8_t parity = raid5_parity_chunk_idx(raid_bdev, stripe_req->stripe_index);
	uint8_t broken = raid_bdev->num_base_bdevs;
	stripe_req_cb cb = NULL;
	int ret = 0;

	for (int i=start; i <= end; i = (i + 1 == parity ? i + 2 : i + 1)) {
		base_ch = raid_bdev_channel_get_base_channel(raid_io->raid_ch, i);
		if (base_ch == NULL) {
			broken = i;
			break;
		}
	}

	if (broken == raid_bdev->num_base_bdevs) {
		// just read (read)
		ret = raid5_map_iovecs_read(stripe_req);
		if (ret) {
			return ret;
		}
		cb = NULL;
	} else {
		// broken req (read_reconstruct)
		ret = raid5_map_iovecs_read_reconstruct(stripe_req, broken);
		if (ret) {
			return ret;
		}
		cb = raid5_read_reconstruct_cb;
	}

	raid5_submit_stripe_request(stripe_req, READ, cb);

	return 0;
}

// rw request submitting start functions

static inline void
raid5_stripe_request_init(
	struct stripe_request *stripe_req,
	struct raid_bdev_io *raid_io,
	uint64_t stripe_index)
{
	stripe_req->raid_io = raid_io;
	stripe_req->stripe_index = stripe_index;

	raid5_chunks_disable(stripe_req);
}

static void
raid5_submit_rw_request(struct raid_bdev_io *raid_io)
{
	struct raid_bdev *raid_bdev = raid_io->raid_bdev;
	struct raid5_info *r5_info = raid_bdev->module_private;
	uint64_t stripe_index = raid_io->offset_blocks / r5_info->stripe_blocks;
	struct raid5_io_channel *r5ch = raid_bdev_channel_get_module_ctx(raid_io->raid_ch);
	struct stripe_request *stripe_req;
	int ret;

	stripe_req = raid5_get_stripe_request(r5ch);
	if (!stripe_req) {
		raid_bdev_io_complete(raid_io, SPDK_BDEV_IO_STATUS_NOMEM);
		return;
	}

	raid5_stripe_request_init(stripe_req, raid_io, stripe_index);

	switch (raid_io->type) {
	case SPDK_BDEV_IO_TYPE_READ:
		ret = raid5_submit_read_request(stripe_req);
		break;
	case SPDK_BDEV_IO_TYPE_WRITE:
		ret = raid5_submit_write_request(stripe_req);
        break;
	default:
		ret = -EINVAL;
		break;
	}

	if (spdk_unlikely(ret)) {
		stripe_req->cb = NULL;
		raid5_stripe_request_complete(stripe_req, raid5_return_code_to_io_status(ret));
	}
}

// raid5 channel creating and destroying function

static void
raid5_stripe_request_free(struct stripe_request *stripe_req)
{
	struct chunk *chunk;

	FOR_EACH_CHUNK(stripe_req, chunk) {
		free(chunk->iovs);
	}

	{
		struct raid5_info *r5_info = raid5_ch_to_r5_info(stripe_req->r5ch);
		struct raid_bdev *raid_bdev = r5_info->raid_bdev;
        uint8_t num_base_bdevs = raid_bdev->num_base_bdevs;
		uint8_t i;

		if (stripe_req->buff.chunk_buffers) {
			for (i = 0; i < num_base_bdevs; i++) {
				if (stripe_req->buff.chunk_buffers[i]) {
					spdk_dma_free(stripe_req->buff.chunk_buffers[i]);
				} else {
					break;
				}
			}
			free(stripe_req->buff.chunk_buffers);
		}

		if (stripe_req->buff.zero_buffer) {
			spdk_dma_free(stripe_req->buff.zero_buffer);
		}
	}

	if (stripe_req->extra.iov_base) {
		spdk_dma_free(stripe_req->extra.iov_base);
	}
    
    free(stripe_req->xor.chunk_iovcnt);
	free(stripe_req->xor.chunk_iovs);
	free(stripe_req->xor.chunk_buffers);
	free(stripe_req->xor.chunk_iov_iters);

	free(stripe_req);
}

static struct stripe_request *
raid5_stripe_request_alloc(struct raid5_io_channel *r5ch)
{
	struct raid5_info *r5_info = raid5_ch_to_r5_info(r5ch);
	struct raid_bdev *raid_bdev = r5_info->raid_bdev;
    uint8_t num_base_bdevs = raid_bdev->num_base_bdevs;
	struct stripe_request *stripe_req;
	struct chunk *chunk;
	size_t buffer_len = raid_bdev->strip_size * raid_bdev->bdev.blocklen;

	stripe_req = calloc(1, sizeof(*stripe_req) + sizeof(*chunk) * raid_bdev->num_base_bdevs);
	if (!stripe_req) {
		return NULL;
	}

	stripe_req->r5ch = r5ch;

	FOR_EACH_CHUNK(stripe_req, chunk) {
		chunk->index = chunk - stripe_req->chunks;
		chunk->iovcnt_max = 4;
		chunk->iovs = calloc(chunk->iovcnt_max, sizeof(chunk->iovs[0]));
		if (!chunk->iovs) {
			goto err;
		}
	}

	{
		void *buf;
		uint8_t i;

		stripe_req->buff.chunk_buffers = calloc(num_base_bdevs, sizeof(void *));
		if (!stripe_req->buff.chunk_buffers) {
			goto err;
		}

		for (i = 0; i < num_base_bdevs; i++) {
			stripe_req->buff.chunk_buffers[i] = NULL;
			buf = spdk_dma_malloc(buffer_len, r5_info->buf_alignment, NULL);
			if (!buf) {
				goto err;
			}
			stripe_req->buff.chunk_buffers[i] = buf;
		}

		stripe_req->buff.zero_buffer = spdk_dma_zmalloc(buffer_len, r5_info->buf_alignment, NULL);
		if (!stripe_req->buff.zero_buffer) {
			goto err;
		}
	}

	stripe_req->extra.iov_base = spdk_dma_malloc(buffer_len, r5_info->buf_alignment, NULL);
	if (!stripe_req->extra.iov_base) {
		goto err;
	}
	stripe_req->extra.iov_len = buffer_len;

    stripe_req->xor.chunk_iov_iters = malloc(SPDK_IOVITER_SIZE(1 + num_base_bdevs));
	if (!stripe_req->xor.chunk_iov_iters) {
		goto err;
	}

    stripe_req->xor.chunk_buffers = calloc(1 + num_base_bdevs, sizeof(*stripe_req->xor.chunk_buffers));
	if (!stripe_req->xor.chunk_buffers) {
		goto err;
	}

	stripe_req->xor.chunk_iovs = calloc(1 + num_base_bdevs, sizeof(*stripe_req->xor.chunk_iovs));
	if (!stripe_req->xor.chunk_iovs) {
		goto err;
	}

	stripe_req->xor.chunk_iovcnt = calloc(1 + num_base_bdevs, sizeof(*stripe_req->xor.chunk_iovcnt));
	if (!stripe_req->xor.chunk_iovcnt) {
		goto err;
	}

	return stripe_req;
err:
	raid5_stripe_request_free(stripe_req);
	return NULL;
}

static void
raid5_ioch_destroy(void *io_device, void *ctx_buf)
{
	struct raid5_io_channel *r5ch = ctx_buf;
	struct stripe_request *stripe_req;

	while ((stripe_req = TAILQ_FIRST(&r5ch->free_stripe_requests))) {
		TAILQ_REMOVE(&r5ch->free_stripe_requests, stripe_req, link);
		raid5_stripe_request_free(stripe_req);
	}
	spdk_spin_destroy(&r5ch->lock);

	if (r5ch->accel_ch) {
		spdk_put_io_channel(r5ch->accel_ch);
	}
}

static int
raid5_ioch_create(void *io_device, void *ctx_buf)
{
	struct raid5_io_channel *r5ch = ctx_buf;
	struct raid5_info *r5_info = io_device;
	struct stripe_request *stripe_req;
	int i;

	TAILQ_INIT(&r5ch->free_stripe_requests);

	for (i = 0; i < RAID5_MAX_STRIPES; i++) {
		stripe_req = raid5_stripe_request_alloc(r5ch);
		if (!stripe_req) {
			goto err;
		}

		TAILQ_INSERT_HEAD(&r5ch->free_stripe_requests, stripe_req, link);
	}
	spdk_spin_init(&r5ch->lock);

	r5ch->accel_ch = spdk_accel_get_io_channel();
	if (!r5ch->accel_ch) {
		SPDK_ERRLOG("Failed to get accel framework's IO channel\n");
		goto err;
	}

	return 0;
err:
	SPDK_ERRLOG("Failed to initialize io channel\n");
	raid5_ioch_destroy(r5_info, r5ch);
	return -ENOMEM;
}

// raid5 starting function

static int
raid5_start(struct raid_bdev *raid_bdev)
{
	uint64_t min_blockcnt = UINT64_MAX;
	uint64_t base_bdev_data_size;
	struct raid_base_bdev_info *base_info;
	struct spdk_bdev *base_bdev;
	struct raid5_info *r5_info;
	size_t alignment = 0;

	r5_info = calloc(1, sizeof(*r5_info));
	if (!r5_info) {
		SPDK_ERRLOG("Failed to allocate r5_info\n");
		return -ENOMEM;
	}
	r5_info->raid_bdev = raid_bdev;

	RAID_FOR_EACH_BASE_BDEV(raid_bdev, base_info) {
		min_blockcnt = spdk_min(min_blockcnt, base_info->data_size);
		if (base_info->desc) {
			base_bdev = spdk_bdev_desc_get_bdev(base_info->desc);
			alignment = spdk_max(alignment, spdk_bdev_get_buf_align(base_bdev));
		}
	}

	base_bdev_data_size = (min_blockcnt / raid_bdev->strip_size) * raid_bdev->strip_size;

	RAID_FOR_EACH_BASE_BDEV(raid_bdev, base_info) {
		base_info->data_size = base_bdev_data_size;
	}

	r5_info->total_stripes = min_blockcnt / raid_bdev->strip_size;
	r5_info->stripe_blocks = raid_bdev->strip_size * raid5_stripe_data_chunks_num(raid_bdev);
	r5_info->buf_alignment = alignment;
	r5_info->stripe_locks = calloc(r5_info->total_stripes, sizeof(struct stripe_lock));
	if(!r5_info->stripe_locks) {
		SPDK_ERRLOG("Failed to allocate stripe_locks\n");
	}
	for (uint64_t i = 0; i < r5_info->total_stripes; i++) {
		spdk_spin_init(&r5_info->stripe_locks[i].lock);
		r5_info->stripe_locks[i].locked = false;
	}

	raid_bdev->bdev.blockcnt = r5_info->stripe_blocks * r5_info->total_stripes;
	raid_bdev->bdev.optimal_io_boundary = r5_info->stripe_blocks;
	raid_bdev->bdev.split_on_optimal_io_boundary = true;

	raid_bdev->module_private = r5_info;

	spdk_io_device_register(r5_info, raid5_ioch_create, raid5_ioch_destroy,
				sizeof(struct raid5_io_channel), NULL);

	return 0;
}

// raid5 stopping funtions

static void
raid5_io_device_unregister_done(void *io_device)
{
	struct raid5_info *r5_info = io_device;
	struct raid_bdev *raid_bdev = r5_info->raid_bdev;

	for (uint64_t i = 0; i < r5_info->total_stripes; i++) {
		spdk_spin_destroy(&r5_info->stripe_locks[i].lock);
	}
	free(r5_info->stripe_locks);

	free(r5_info);

	raid_bdev_module_stop_done(raid_bdev);
}

static bool
raid5_stop(struct raid_bdev *raid_bdev)
{
	struct raid5_info *r5_info = raid_bdev->module_private;

	spdk_io_device_unregister(r5_info, raid5_io_device_unregister_done);

	return false;
}

// raid5 channel acquisition function

static struct spdk_io_channel *
raid5_get_io_channel(struct raid_bdev *raid_bdev)
{
	struct raid5_info *r5_info = raid_bdev->module_private;

	return spdk_get_io_channel(r5_info);
}

// initializing and registering a raid module

static struct raid_bdev_module g_raid5_module = {
	.level = RAID5,
	.base_bdevs_min = 3,
	.base_bdevs_constraint = {CONSTRAINT_MAX_BASE_BDEVS_REMOVED, 1},
	.start = raid5_start,
	.stop = raid5_stop,
	.submit_rw_request = raid5_submit_rw_request,
	.get_io_channel = raid5_get_io_channel,
};
RAID_MODULE_REGISTER(&g_raid5_module)

SPDK_LOG_REGISTER_COMPONENT(bdev_raid5)
