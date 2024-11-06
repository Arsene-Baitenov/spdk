#!/bin/bash

# before running:
# load the module ublk_drv.ko into the kernel: insmod {...}/ublk_drv.ko
# activate root mode: sudo -i
# activate virtualenv: source {virtualenv dir}/bin/activate
# run scipt from spdk dir: ./configure --with-ublk
# type 'make' to build spdk
# run this script: bash {path to script}/run_test.sh {full path to spdk}

# $1 full path to spdk

# exit codes:
# 1 - path to spdk isn't entered
# 2 - ublk target isn't created

spdk=$1;

function start() {
    local old_dir=$(pwd);
    cd $spdk;

    ./scripts/setup.sh >/dev/null&
    pid=$!;
    wait "$pid";

    start-stop-daemon -Sbv -n spdk_tgt -x $spdk/build/bin/spdk_tgt;
    sleep 5;
    ./scripts/rpc.py ublk_create_target;
    if [[ "$?" != "0" ]]; then
        start-stop-daemon -Kvx $spdk/build/bin/spdk_tgt;
        cd $old_dir;
        exit 2;
    fi

    cd $old_dir;
}

function finish() {
    local old_dir=$(pwd);
    cd $spdk;

    ./scripts/rpc.py ublk_destroy_target;
    start-stop-daemon -Kvx $spdk/build/bin/spdk_tgt;
    cd $old_dir;
}

# $1 full path to fio config
function run_test_with_fio() {
    local old_dir=$(pwd);
    local fio_conf=$1;
    local raid_start_conf=$spdk/test/bdev/raid5/spdk-conf/start.json;
    local raid_crash_conf=$spdk/test/bdev/raid5/spdk-conf/crash_base_bdev.json;
    local raid_stop_conf=$spdk/test/bdev/raid5/spdk-conf/stop.json;
    
    cd $spdk;
    ./scripts/rpc.py load_config -j $raid_start_conf;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 1;
    fi

    ./scripts/rpc.py ublk_start_disk Raid5 1 >/dev/null;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 1;
    fi

    fio $fio_conf >/dev/null;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 2;
    fi
    
    ./scripts/rpc.py load_config -j $raid_crash_conf;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 1;
    fi

    fio $fio_conf >/dev/null;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 2;
    fi
    
    ./scripts/rpc.py ublk_stop_disk 1;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 1;
    fi
    
    ./scripts/rpc.py load_config -j $raid_stop_conf;
    if [[ "$?" != "0" ]]; then
        cd $old_dir;
        return 1;
    fi

    cd $old_dir;
    return 0;
}

# $1 prestate
# $2 full path to fio config
# $3 name
function run_test() {
    local prestate=$1;
    local fio_conf=$2;
    local name=$3;

    if [[ "$prestate" != "0" ]]; then
        echo -e "$name: \033[33mskipped\033[0m";
        return $prestate;
    else
        run_test_with_fio $fio_conf;
        if [[ "$?" == 0 ]]; then
            echo -e "$name: \033[32mpassed\033[0m";
            return 0;
        else
            echo -e "$name: \033[31mfailed\033[0m";
            return 1;
        fi
    fi
}

function run_tests() {
    local fio_conf_dir=$spdk/test/bdev/raid5/fio-conf;
    
    run_test "0" $fio_conf_dir/write-small.fio "write small";
    run_test $? $fio_conf_dir/write-big.fio "write big";
    run_test $? $fio_conf_dir/randwrite-small.fio "randwrite small";
    run_test $? $fio_conf_dir/randwrite-big.fio "randwrite big";
    run_test $? $fio_conf_dir/iodepth4.fio "iodepth=4";
    run_test $? $fio_conf_dir/iodepth16.fio "iodepth=16";
}

if [[ -z "$spdk" ]]; then
    echo "error: path to spdk isn't entered"
    exit 1
fi

start;
echo;

run_tests;

echo;
finish;