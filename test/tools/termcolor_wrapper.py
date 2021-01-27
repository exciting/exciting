try:
    from termcolor import colored
    color = True
except ModuleNotFoundError:
    color = False

def print_color(message:str, text_color:str, end='\n'):
    """                                        
    Wrapper to print function, allowing colored text    
    
    :param message: string to print. 
    :param color: string indicating text color.
    :param end: optional string appended after the last value, default a newline.
    """
    if color:
        print(colored(message, text_color), end=end)
    else:
        print(message, end=end)
    return
