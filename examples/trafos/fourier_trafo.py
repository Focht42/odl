# Copyright 2014-2016 The ODL development group
#
# This file is part of ODL.
#
# ODL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ODL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ODL.  If not, see <http://www.gnu.org/licenses/>.

"""Simple example on the usage of the Fourier Transform."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

# Internal
import odl


# Discretized space: discretized functions on the rectangle [-1, 1] x [-1, 1]
# with 512 samples per dimension and complex data type (for full FT).
space = odl.uniform_discr([-1, -1], [1, 1], (512, 512), dtype='complex')

# Make the Fourier transform operator on this space. The range is calculated
# automatically. The default backend is numpy.fft
ft_op = odl.trafos.FourierTransform(space)

# Create a phantom and its Fourier transfrom and display them
phantom = odl.phantom.shepp_logan(space, modified=True)
phantom.show(title='Shepp-Logan phantom')
phantom_ft = ft_op(phantom)
phantom_ft.show(title='Full Fourier Transform')

# Calculate the inverse transform
phantom_ft_inv = ft_op.inverse(phantom_ft)
phantom_ft_inv.show(title='Full Fourier Transform inverted')

# Calculate the FT only along the first axis
ft_op_axis0 = odl.trafos.FourierTransform(space, axes=0)
phantom_ft_axis0 = ft_op_axis0(phantom)
phantom_ft_axis0.show(title='Fourier transform along axis 0')

# If a real space is used, the Fourier transform can be calculated in the
# "half-complex" mode. This means that along the last axis of the transform,
# only the negative half of the spectrum is stored since the other half is
# its complex conjugate. This is faster and more memory efficient.
real_space = space.real_space
ft_op_halfc = odl.trafos.FourierTransform(real_space, halfcomplex=True)
phantom_real = odl.phantom.shepp_logan(real_space, modified=True)
phantom_real.show(title='Shepp-Logan phantom, real version')
phantom_real_ft = ft_op_halfc(phantom_real)
phantom_real_ft.show(title='Half-complex Fourier Transform')

# If the space is real, the inverse also gives a real result
phantom_real_ft_inv = ft_op_halfc.inverse(phantom_real_ft)
phantom_real_ft_inv.show(title='Half-complex Fourier Transform inverted',
                         show=True)
