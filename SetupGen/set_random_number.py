import os.path
import random


def log_rand(rand):
    logger_path = "rng_logger.txt"
    if os.path.exists(logger_path):
        with open(logger_path, "a") as f:
            f.write(f"set new randon num {rand} \n")
        f.close()
        return
    with open(logger_path, "w") as f:
        f.write(f"set new randon num {rand} \n")
        f.close()


with open("gpu_step5_production.inp", "r") as f:
    lines = f.readlines()
f.close()

for i, line in enumerate(lines):
    if "set SVALUE" in line:
        rv = random.randint(147483648, 2100083641)
        lines[i] = f"set SVALUE {rv}\n"
        log_rand(rv)
    if "set equi 1" in line:
        lines[i] = "set equi 0\n"
        log_rand("EQUI")

with open("gpu_step5_production.inp", "w") as f:
    f.writelines(lines)
f.close()
