def copy_path(root_src_dir, root_dst_dir, exclude_files=[]):

    from shutil import copy
    import os

    for src_dir, dirs, files in os.walk(root_src_dir):
        dst_dir = src_dir.replace(root_src_dir, root_dst_dir)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for f in files:
            if f not in exclude_files:
                src_file = os.path.join(src_dir, f)
                dst_file = os.path.join(dst_dir, f)
                copy(src_file, dst_file)
