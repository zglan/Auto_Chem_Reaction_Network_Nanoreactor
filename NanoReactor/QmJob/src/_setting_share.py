from ._setting_default import default_setting

class SharedSetting:
    def __init__(self):
        for key in default_setting:
            setattr(self, key, default_setting[key])
