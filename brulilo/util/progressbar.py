"""
Simple class for a text-based progress bar.  Nothing too fancy here
"""
import sys


class ProgressBar(object):
    _bar_width = 60
    _p_symbol = "#"
    _r_symbol = "-"
    _done_status = "Done!\n"

    def __init__(self, label):
        self.label = label
        self.state = 0
        self._update_display()

    def update(self, status=''):
        self.state += self.delta_update
        self._update_display(status)

    def _update_display(self, status=''):
        progress = int(self._bar_width*self.state)
        remaining = self._bar_width - progress
        if self.state == 1.0:
            status = self._done_status
        display = "\r{0}: [{1}] {2}".format(self.label,
                                            (self._p_symbol*progress +
                                             self._r_symbol*remaining),
                                            status)
        sys.stdout.write(display)
        sys.stdout.flush()


class IntProgressBar(ProgressBar):
    def __init__(self, label='Progress', num_ticks=1):
        self.total_update = num_ticks
        self.delta_update = 1./self.total_update

        super(IntProgressBar, self).__init__(label)
