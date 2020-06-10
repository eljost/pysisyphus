class NeedNewInternalsException(Exception):
    def __init__(self, cart_coords, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.cart_coords = cart_coords


class RebuiltInternalsException(Exception):
    pass
