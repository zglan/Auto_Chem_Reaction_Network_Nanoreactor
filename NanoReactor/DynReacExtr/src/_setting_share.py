from _setting_default import default_setting
from _setting_arg import ArgParse


class SharedSetting:
    def __init__(self):
        args = vars(ArgParse.get_args().parse_args())
        for key in args:
            setattr(self, key, args[key])
        for key in default_setting:
            setattr(self, key, default_setting[key])