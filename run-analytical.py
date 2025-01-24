import os, argparse

# Option parser
parser = argparse.ArgumentParser(description="MC/DC Verification - Analytical")
parser.add_argument("--srun", type=int, default=0)
parser.add_argument("--mpiexec", type=int, default=0)
parser.add_argument("--mpirun", type=int, default=0)
args, unargs = parser.parse_known_args()

# Get the MPI option
if args.srun > 0 or args.mpiexec > 0 or args.mpirun > 0:
    if args.srun > 1:
        mpi_option = "--srun=%i" % args.srun
    elif args.mpiexec > 1:
        mpi_option = "--mpiexec=%i" % args.mpiexec
    elif args.mpirun > 1:
        mpi_option = "--mpirun=%i" % args.mpirun

# Analytical verification
os.chdir("analytical")
os.system("python run.py %s" % mpi_option)
