open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    printf '%.3d' $? >&3
    )&
}

N=13
open_sem $N
for i in {0..127}
do
	for j in {5..7}
	do
    	run_with_lock ./Space4ErrBody primaryTHESIS_HORUS.txt THESISnodalParameters_${j}nodes.txt Calibration ${i} > log_THESIS_HORUS_${j}nodes_Calibration_${i}.dat 2>&1 & disown;
	done 
done