"""
Module defining custom wildcards, and their corresponding regex replacements
This redefines regex pattern matches with more compact notation.
"""
from collections import namedtuple

Wildcards = namedtuple('Wildcard', ['word',
                                    'character',
                                    'number',
                                    'letter'])

wildcard = Wildcards(word='??',  # wildcard for string of letters, digits and '-'
                     character='*',  # wildcard for single letter or digit
                     number='#',  # wildcard for digit
                     letter='&'  # wildcard for letter
                     )

regex = Wildcards(word=r'[\w-]+',  # regular expression for string of letters, digits and '-'
                  character=r'[a-zA-Z\d]',  # regular expression for letter or digit
                  number=r'[\d]',  # regular expression for digit
                  letter=r'[a-zA-Z]'  # regular expression for letter
                  )


def wildcard_processor(string: str) -> str:
    """
    Replace wildcards in a string with regular expressions as defined by wildcard and regex,
    replace '.' with '[.]' and return the regular expression pattern.
    Examples:
      'abc.def' --> '^abc[.]def'
      'abc_??' --> '^abc_[\w-]$'
      'abc_*'  --> '^abc_[a-zA-Z\d]$'
      'abc_#'  --> '^abc_[\d]$'
      'abc_&'  --> '^abc_[a-zA-Z]$'
      '*abc_??_sdf*_##_qwe.oiu_&.out'
               --> '^[a-zA-Z\d]$abc_[\w-]_sdf[a-zA-Z\d]_[\d][\d]_qwe[.]oiu_[a-zA-Z][.]out$'

    :param string: string with wildcards
    :return str: regular expressions pattern
    """
    re_string = '^' + string + '$'
    re_string = re_string.replace('.', '[.]')
    re_string = re_string.replace(wildcard.word, regex.word)
    re_string = re_string.replace(wildcard.character, regex.character)
    re_string = re_string.replace(wildcard.number, regex.number)
    re_string = re_string.replace(wildcard.letter, regex.letter)

    return re_string
