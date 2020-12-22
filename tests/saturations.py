import subprocess
import glob
from threading import Thread

files = glob.glob('settings/*')
files = [f for f in files if 'um' in f]

def run(fname):
    subprocess.run(["python","saturation.py",f"-f{fname}","-D","-NDCR","-NXT","-NAP"])

threads = []
for f in files:
    t = Thread(target=run, args=(f,))
    t.start()
    threads.append(t)

for t in threads:
    t.join()
