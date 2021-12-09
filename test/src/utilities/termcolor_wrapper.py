try:
    from termcolor import colored

    def print_color(message: str, text_color: str, end='\n'):
        """
        Wrapper to print function, allowing colored text

        :param message: string to print.
        :param text_color: string indicating text color.
        :param end: optional string appended after the last value, default a newline.
        """
        print(colored(message, text_color), end=end)

except ModuleNotFoundError:
    def print_color(message: str, text_color: str, end='\n'):
        print(message, end=end)
