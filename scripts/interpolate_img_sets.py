'''
Interpolates sets of images from h5 files to higher resolution for publication-quality
images
Usage:
    interpolate_img_sets.py <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames]
'''
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.ioff()
from dedalus.extras import plot_tools

def main(filename, start, count, output, _min, _max):
    dpi = 200
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_{:06}.png'.format(write)
    tasks = ['density']

    with h5py.File(filename, mode='r') as file:
        fig = plt.figure()
        for index in range(start, start+count):
            for n, task in enumerate(tasks):
                dset = np.array(file['tasks'][task][index])

                psm = plt.imshow(np.rot90(dset), cmap='jet', interpolation='bicubic', vmin=_min, vmax=_max)
                plt.title(title_func(file['scales/sim_time'][index]))

                f2 = plt.gca()
                f2.axes.get_xaxis().set_visible(False)
                f2.axes.get_yaxis().set_visible(False)

            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
        plt.close(fig)

if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync
    import glob

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()

    post.visit_writes(args['<files>'], main, output=output_path, _min=-0.5, _max=0.5)
