import os
import shutil


def clear_files(event_path):
    files = ['boat_fault.txt', 'dyfi_dat.xml', 'event.xml', 'model.conf',
             'model_select.conf', 'moment.xml', 'stationlist.xml',
             'shake_data.hdf', 'dir_fault.txt']
    for fname in files:
        try:
            os.remove(os.path.join(event_path, fname))
        except OSError:
            pass


def set_files(event_path, files):
    for src, dst in files.items():
        shutil.copy(os.path.join(event_path, 'data', src),
                    os.path.join(event_path, dst))
