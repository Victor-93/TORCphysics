
# This class will be use to create different type of events
class Event:
    def __init__(self, time, frame, event_type, twist, superhelical, global_twist, global_superhelical,
                 site, enzyme, position):
        self.time = time
        self.frame = frame
        self.event_type = event_type
        self.twist = twist
        self.superhelical = superhelical
        self.global_twist = global_twist
        self.global_superhelical = global_superhelical
        self.site = site
        self.enzyme = enzyme
        self.position = position

