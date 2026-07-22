def list_categories_and_scripts(path_):
    import os
    dirname = os.path.dirname(path_)
    dir_content = os.listdir(dirname)
    scripts = []
    categories = []
    for file_ in sorted(dir_content):
        if os.path.isdir(os.path.join(dirname, file_)) and not file_.startswith('__'):
            if "__main__.py" in os.listdir(os.path.join(dirname, file_)):
                categories.append(os.path.basename(file_))
                continue
        elif file_.endswith('.py') and not file_.startswith('_'):
            scripts.append(os.path.basename(file_).removesuffix('.py'))
    if len(categories)>0:
        print('Available categories:')
        for cat in categories:
                print(f"- {cat}")
        print()
    print('Available scripts in this category:')
    for script in scripts:
            print(f"- {script}")