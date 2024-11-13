ifort global_parameters.f90 utility_mod.f90 main.f90 -o test
for f in 0
# for f in {1..10}
do
    mkdir /home/zhangjing/q0/$f
	cd /home/zhangjing/q0/$f
	cp /home/zhangjing/q0/input.txt /home/zhangjing/q0/$f
	cp /home/zhangjing/q0/tql.pbs /home/zhangjing/q0/$f
	cp /home/zhangjing/q0/test /home/zhangjing/q0/$f
	cp /home/zhangjing/q0/delta.txt /home/zhangjing/q0/$f
	echo -e "$f" >> /home/zhangjing/q0/"$f"/input.txt

# 	# cat<<EOF>input.txt
# # $f
# # EOF
	qsub tql.pbs
    cd ..
done
