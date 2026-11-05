import types
import unittest
from unittest import mock

import numpy as np
import xarray as xr

from src import preprocessor


def _var(name, lat=np.arange(3.0)):
    # VarlistEntry stand-in: translation is on the model's yh/xh, the POD wants lat/lon
    dim_axes = {'Y': types.SimpleNamespace(axis='Y', name='yh'),
                'X': types.SimpleNamespace(axis='X', name='xh')}
    return types.SimpleNamespace(
        name=name, full_name=f'<case.{name}>', log=mock.MagicMock(),
        translation=types.SimpleNamespace(name=name, dim_axes=dim_axes, scalar_coords=[]),
        axes={'Y': types.SimpleNamespace(name='lat'), 'X': types.SimpleNamespace(name='lon')}
    )


def _dataset(yh=np.arange(3.0)):
    # one shared dataset, as in process(): 'so' was already renamed onto lat/lon,
    # 'thetao' was re-added from the unrenamed data and is still on yh/xh
    lon = np.arange(4.0)
    return xr.Dataset(
        {'so': (('lat', 'lon'), np.ones((3, 4))),
         'thetao': (('yh', 'xh'), np.zeros((len(yh), 4)))},
        coords={'lat': np.arange(3.0), 'lon': lon, 'yh': yh, 'xh': lon}
    )


class TestRenameOntoExistingCoords(unittest.TestCase):
    def _execute(self, var, ds):
        func = object.__new__(preprocessor.RenameVariablesFunction)
        return func.execute(var, ds)

    def test_second_variable_moved_onto_renamed_coords(self):
        ds = self._execute(_var('thetao'), _dataset())
        self.assertEqual(ds['thetao'].dims, ('lat', 'lon'))
        self.assertEqual(ds['so'].dims, ('lat', 'lon'))
        self.assertNotIn('yh', ds.coords)
        self.assertNotIn('xh', ds.coords)

    def test_coords_with_different_values_left_alone(self):
        var = _var('thetao')
        ds = self._execute(var, _dataset(yh=np.arange(3.0) + 0.5))
        self.assertEqual(ds['thetao'].dims, ('yh', 'lon'))
        self.assertIn('yh', ds.coords)
        var.log.warning.assert_called_once()

    def test_first_variable_renamed_as_before(self):
        ds = xr.Dataset({'thetao': (('yh', 'xh'), np.zeros((3, 4)))},
                        coords={'yh': np.arange(3.0), 'xh': np.arange(4.0)})
        ds = self._execute(_var('thetao'), ds)
        self.assertEqual(ds['thetao'].dims, ('lat', 'lon'))


if __name__ == '__main__':
    unittest.main()
