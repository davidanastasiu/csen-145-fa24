EXECS=hello_world hello_world_omp check_status probe groups split my_bcast compare_bcast reduce_avg reduce_stddev avg all_avg send_recv ping_pong ring bin random_rank random_walk
MPICC?=mpicc
MPICXX?=mpicxx

all: ${EXECS}

hello_world: hello_world.c
	${MPICC} -o hello_world hello_world.c

hello_world_omp: hello_world_omp.c
	${MPICC} -fopenmp -o hello_world_omp hello_world_omp.c

check_status: check_status.c
	${MPICC} -o check_status check_status.c

probe: probe.c
	${MPICC} -o probe probe.c

split: split.c
	${MPICC} -o split split.c

groups: groups.c
	${MPICC} -o groups groups.c

my_bcast: my_bcast.c
	${MPICC} -o my_bcast my_bcast.c

compare_bcast: compare_bcast.c
	${MPICC} -o compare_bcast compare_bcast.c

reduce_avg: reduce_avg.c
	${MPICC} -o reduce_avg reduce_avg.c

reduce_stddev: reduce_stddev.c
	${MPICC} -o reduce_stddev reduce_stddev.c -lm

avg: avg.c
	${MPICC} -o avg avg.c

all_avg: all_avg.c
	${MPICC} -o all_avg all_avg.c

send_recv: send_recv.c
	${MPICC} -o send_recv send_recv.c

ping_pong: ping_pong.c
	${MPICC} -o ping_pong ping_pong.c

ring: ring.c
	${MPICC} -o ring ring.c

bin: bin.c
	${MPICC} -o bin bin.c

tmpi_rank.o: tmpi_rank.c
	${MPICC} -c tmpi_rank.c

random_rank: tmpi_rank.o random_rank.c
	${MPICC} -o random_rank random_rank.c tmpi_rank.o

random_walk: random_walk.c
	${MPICXX} -o random_walk random_walk.c

clean:
	rm -f ${EXECS}
