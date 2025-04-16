import os
import shutil


def get_subfolders(path):
    folder_paths = []
    if os.path.isdir(path):
        _folders = os.listdir(path)
        for _folder in _folders:
            _fpath = os.path.join(path, _folder)
            if os.path.isdir(_fpath):
                folder_paths.append(_fpath)

    return folder_paths


cwd = os.getcwd()
folders = os.listdir(cwd)
new_observer_path = os.path.join(cwd, "Observer.py")


depth = 2


if depth == 1:
    for folder in folders:
        if os.path.isdir(folder):
            folder_path = os.path.join(cwd, folder)
            observer_path = os.path.join(folder_path, "Observer.py")
            if os.path.exists(observer_path):
                os.remove(observer_path)
                shutil.copy(new_observer_path, observer_path)
                print("replacing", observer_path)

if depth == 2:

    print("starting replacing at recursive depth 2")

    for folder in folders:
        folder_list = get_subfolders(os.path.join(os.getcwd(), folder))
        for f_path in folder_list:
            observer_path = os.path.join(f_path, "Observer.py")
            if os.path.exists(observer_path):
                os.remove(observer_path)
                shutil.copy(new_observer_path, observer_path)
                print("replacing", observer_path)