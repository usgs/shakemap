#!/usr/bin/env python
import re
import warnings
from datetime import datetime


class Specifier(str):
    """Model %Y and such in `strftime`'s format string."""
    def __new__(cls, *args):
        self = super(Specifier, cls).__new__(cls, *args)
        assert self.startswith('%')
        assert len(self) == 2
        self._regex = re.compile(r'(%*{0})'.format(str(self)))
        return self

    def ispresent_in(self, format):
        m = self._regex.search(format)
        return m and m.group(1).count('%') & 1  # odd number of '%'

    def replace_in(self, format, by):
        def repl(m):
            n = m.group(1).count('%')
            if n & 1:  # odd number of '%'
                prefix = '%' * (n - 1) if n > 0 else ''
                return prefix + str(by)  # replace format
            else:
                return m.group(0)  # leave unchanged
        return self._regex.sub(repl, format)


class ShakeDateTime(datetime):

    def strftime(self, format):
        year = self.year
        if year >= 1900:
            return super(ShakeDateTime, self).strftime(format)
        assert year < 1900
        factor = (1900 - year - 1) // 400 + 1
        future_year = year + factor * 400
        assert future_year > 1900

        format = Specifier('%Y').replace_in(format, year)
        result = self.replace(year=future_year).strftime(format)
        if any(f.ispresent_in(format) for f in map(Specifier, ['%c', '%x'])):
            msg = "'%c', '%x' produce unreliable results for year < 1900"
            if not force:
                raise ValueError(msg + " use force=True to override")
            warnings.warn(msg)
            result = result.replace(str(future_year), str(year))
        assert (future_year % 100) == (year %
                                       100)  # last two digits are the same
        return result

if __name__ == "__main__":
    t1 = ShakeDateTime(1900, 1, 1, 1, 1, 1)
    print(t1.strftime('%Y-%m-%d %H:%M:%S'))
    t2 = ShakeDateTime(1899, 7, 13, 12, 16, 41)
    print(t2.strftime('%Y-%m-%d %H:%M:%S'))
    t3 = ShakeDateTime.strptime('1899-01-01', '%Y-%m-%d')
    print(t3.strftime('%Y-%m-%d'))
