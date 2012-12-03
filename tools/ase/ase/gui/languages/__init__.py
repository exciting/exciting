import locale

from ase.gui.languages.en import translation as default_translation


def translate(text):
    return default_translation.get(text, text)

language_code = locale.getdefaultlocale()[0]
if language_code is None:
    language_code = 'en'
else:
    language_code = language_code[:2]

if language_code != 'en':
    try:
        module = __import__(language_code, globals(), locals())
    except ImportError:
        pass
    else:
        translation = module.translation
        def translate(text):
            translated_text = translation.get(text)
            if translated_text is None:
                return default_translation.get(text, text)
            else:
                return translated_text
