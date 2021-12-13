import sys
import subprocess
from collections import namedtuple
from matplotlib import pyplot as plt

FOMS = namedtuple('foms', 'crcs, snrs')
CRC  = namedtuple('crc', 'crc, var')

iterations = (0, 4, 9, 14, 19, 24, 29)
iterations = (0, 14, 29)
#iterations = tuple(range(30))

xxx = dict(true='Truelor', reco='mlor')
TRUE, RECO = 'Truemlor', 'mlor'

def tof_name(tof): return 'tof{:->3}'.format(tof or '')
def   q_name(q  ): return   'q{:->4}'.format(q   or '')
def   E_name(E  ): return   'E{:->3}'.format(E   or '')

def generate_dir(lors, mm, tof, q, E):
    return f'data/jaszczak/imgs/udlx-s0-r0-b1-vacuum-nosteel-dz1m-noQrz-6/{lors}/LXe{mm}mm/{tof_name(tof)}-{q_name(q)}-{E_name(E)}'

runs = ((TRUE, 20, None, None,  434), #  0
        (TRUE, 20,   49, None,  434), #  1
        (TRUE, 20,   98, None,  434), #  2
        (TRUE, 20, None, 1900, None), #  3
        (TRUE, 20,   49, 1900, None), #  4
        (TRUE, 20,   98, 1900, None), #  5

        (TRUE, 30, None, None,  434), #  6
        (TRUE, 30,   62, None,  434), #  7
        (TRUE, 30,  125, None,  434), #  8
        (TRUE, 30, None, 1900, None), #  9
        (TRUE, 30,   62, 1900, None), # 10
        (TRUE, 30,  125, 1900, None), # 11

        (TRUE, 40, None, None,  434), # 12
        (TRUE, 40,   79, None,  434), # 13
        (TRUE, 40,  157, None,  434), # 14
        (TRUE, 40, None, 1700, None), # 15
        (TRUE, 40,   79, 1700, None), # 16
        (TRUE, 40,  157, 1700, None), # 17

        (RECO, 20, None, 1900, None), # 18
        (RECO, 20,   49, 1900, None), # 19
        (RECO, 20,   98, 1900, None), # 20

        (RECO, 30, None, 1900, None), # 21
        (RECO, 30,   62, 1900, None), # 22
        (RECO, 30,  125, 1900, None), # 23

        (RECO, 40, None, 1700, None), # 24
        (RECO, 40,   79, 1700, None), # 25
        (RECO, 40,  157, 1700, None)) # 26

sphere_diameters = 9.5, 12.7, 15.9, 19.1, 25.4, 31.8

def get_foms(lors, mm, tof, q, E, it):

    directory = generate_dir(lors, mm, tof, q, E)
    command = f'cargo run --release --bin foms -- {directory}/{it:02}.raw jaszczak'
    o = subprocess.run(command, shell=True, capture_output=True)
    crcs = {}
    snrs = {}

    for line in o.stdout.decode('utf-8').split('\n'):

        columns = line.split()
        if len(columns) == 4:
            d, crc, var, snr = map(float, columns)
            crcs[d] = CRC(crc=crc, var=var)
            snrs[d] = snr

    return FOMS(crcs=crcs, snrs=snrs)


from contextlib import contextmanager
@contextmanager
def dummy(*args, **kwds):
    yield

for run_spec in runs[23:24]:

    def write(*args, **kwds):
        print(*args, **kwds, file=sys.stdout)
        print(*args, **kwds, file=outfile)

    directory = generate_dir(*run_spec)
    filename = f'{directory}/foms'
    print(filename)
    with open(filename, 'w') as outfile:
        write('  ', end='')
        ds_header = ''.join(f'{d:7.1f}' for d in sphere_diameters)
        write(f'{"CRCs":>25}        {"background variabilities":>53}          {"SNRs":>31}', end='\n\n')
        write(f'  {ds_header}         {ds_header}         {ds_header}')
        write()
        for it in range(30):
            img_spec = run_spec + (it,)
            foms = get_foms(*img_spec)
            crcs, snrs = foms
            write(f'{it+1:2} ', end='')
            crcs, variabilities = tuple(zip(*((crc, var) for (r, (crc, var)) in crcs.items())))
            write(''.join(f'{c:6.1f} ' for c in crcs)         , end='         ')
            write(''.join(f'{v:6.1f} ' for v in variabilities), end='         ')
            write(''.join(f'{r:6.1f} ' for r in snrs.values())) # look broken in the Rust implementation of foms


#TODO add ability to read foms from generated fom files

exit(0)

fs = {}

run_numbers_to_use = (20,)

for run_number in run_numbers_to_use:
    for it in iterations:
        run_spec = runs[run_number]
        img_spec = run_spec + (it,)
        print(img_spec)
        print(get_foms(*img_spec))
        fs[img_spec] = get_foms(*img_spec)

for d in sphere_diameters:
    plt.figure()
    count = 0
    for run_number in run_numbers_to_use:
        run_spec = runs[run_number]
        img_spec = run_spec + (it,)
        y = [fs[img_spec].crcs[d].crc for i in iterations]
        e = [fs[img_spec].crcs[d].var for i in iterations]
        x = [1 + n - count/10 for n in iterations]
        #x = iterations
        #plt.errorbar(x, y, yerr=e, label=f'TOF={t} {c}', capsize=3)
        plt.errorbar(x, y, yerr=e, label=f'TOF=t c', capsize=3)
        count += 1
    plt.legend()
    plt.ylim(top=120)
    plt.title(f'CRC vs iteration for {d}mm sphere')
    plt.savefig(f'crcs-{d}mm.png')

plt.show()
