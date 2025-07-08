import importlib
import inspect
import pydoc
import pkgutil
import os
import re

import excitingscripts

TEMPLATE_NAME = "README.template"

def docstring_to_markdown(docstring: str) -> str:
    '''Replace newlines with two spaces and a newline, to avoid wrapping of lines in markdown.'''
    return re.sub('\n', '  \n', docstring)

def extract_docs():
    output = ""
    # iterate through modules of excitingscripts
    for module_info in pkgutil.walk_packages(excitingscripts.__path__,  excitingscripts.__name__ + '.'):
        # if the module is a package, generate a header and print docs if they exist
        if module_info.ispkg:
            output += f"### {module_info.name}\n\n"
            module_docs = pydoc.getdoc(module_info.name)
            if len(module_docs) > 0:
                output += f"{docstring_to_markdown(module_docs)}\n\n"
        # the module is a file with code
        else:
            output += f"#### {module_info.name}\n\n"
            # import the module and iterate over its members (i.e. functions, classes, or data)
            module = importlib.import_module(module_info.name)
            module_docs = pydoc.getdoc(module)
            if len(module_docs) > 0:
                output += f"{docstring_to_markdown(module_docs)}\n\n"
            for mem in inspect.getmembers(module):
                # ignore the main (because it is only used to indicate what is executed when the module is called as a script)
                if mem[0] == 'main':
                    continue
                # if the members path is from this module of excitingscripts, extract the docs and print them here
                try:
                    module_path = inspect.getmodule(mem[1]).__name__.split('.')
                    if module_path[0] == 'excitingscripts' and module_path[:-1] == module_info.name.split('.')[:-1]:
                        output += f"##### {mem[0]}\n\n{docstring_to_markdown(pydoc.getdoc(mem[1]))}\n\n"
                except:  # noqa: E722
                    pass
        output += "\n"
    return output

def main():
    with open(TEMPLATE_NAME, 'r') as f_:
        template = f_.read()

    readme_text = template + extract_docs()

    if not os.path.exists('README.md') or True:
        with open('README.md', 'w') as f_:
            f_.write(readme_text)
    else:
        raise FileExistsError("README.md is already present in this directory. Please delete or rename.")

if __name__ == '__main__':
    
    main()