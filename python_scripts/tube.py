import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from typing import List

"""
# Definition of tube class
"""

class fiber:
  """
  class describing a cnt fiber
  """

  def __init__(self, r, chiral, scaleFactor=10):
    '''
    Class constructor
    '''
    self._r = np.array(r)
    self._scaleFactor = scaleFactor
    self.scale(scaleFactor)
    self.chiral = chiral

  def get_chiral(self):
    return self.chiral

  def r(self, mode='rough', n=100, s=500, k=3):
    '''
    get coordinates of fiber in 'rough' or 'fine' mode

    Parameters:
      mode (str): determines if 'rough' mesh is returned or 'fine' mesh is returned
      n (int): integer indicating the number of points in the fine mesh
      s (int): A smoothing condition. Larger s means more smoothing while smaller values of s indicate less smoothing.
      k (int): Degree of spline
    '''
    if mode == 'rough':
      return self._r
    if mode == 'fine':
      if not hasattr(self, '_r_fine'):
        self._r_fine = self.calculate_r_fine(n=n, s=s, k=k)
      return self._r_fine
    raise 'Unknown mode!'

  def num_nodes(self, mode='rough'):
    '''
    Get the number of nodes in the rough mesh
    '''
    return self.r(mode=mode).shape[0]

  def calculate_r_fine(self, n=100, s=500, k=3):
    """
    Interpolate the rough curve described by self._r to get finer mesh points
    interpolation uses B-spline curves
    
    Parameters:
      n (int): integer indicating the number of points in the fine mesh
      s (int): A smoothing condition. Larger s means more smoothing while smaller values of s indicate less smoothing.
      k (int): Degree of spline

    Return:
      np.ndarray of shape (n,3) where n is the new number of points
    """
    tck, u = interpolate.splprep([self._r[:,0], self._r[:,1], self._r[:,2]], s=s, k=3)
    u_fine = np.linspace(0,1,n)
    x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)
    r_fine = np.stack((x_fine, y_fine, z_fine), axis=-1)
    self._tck, self._u = tck, u
    return r_fine

  def r_fine(self, n=100, s=500, k=3):
    '''
    get refined mesh of fiber coordinates

    Parameters:
      n (int): integer indicating the number of points in the fine mesh
      s (int): A smoothing condition. Larger s means more smoothing while smaller values of s indicate less smoothing.
      k (int): Degree of spline

    Returns:
      np.ndarray of shape (n,3) where n is the number of points in the refined mesh in each fiber
    '''
    if not hasattr(self, '_r_fine'):
      self._r_fine = self.calculate_r_fine(n=n, s=s, k=k)
    return self._r_fine

  def tangent_vec(self) -> np.ndarray:
    """
    Returns:
      np.ndarray of shape (N,3): normalized tangent vectors at the interpolated points
    """
    if not hasattr(self, '_tck'):
      self.calculate_r_fine()
    tck = self._tck
    u_fine = np.linspace(0, 1, self.r(mode='fine').shape[0])
    deriv_x, deriv_y, deriv_z = interpolate.splev(u_fine, tck, der=1)
    t_vec = np.stack((deriv_x, deriv_y, deriv_z), axis=1)
    
    normalize = lambda v: v/np.linalg.norm(v)
    t_vec = np.apply_along_axis(normalize, axis=1, arr=t_vec)
    return t_vec

  def scale(self, factor=10):
    '''
    Scale the coordinates of the rough fiber mesh by a multiplication factor

    Note:
      We do not change the orientation of the sections
    '''
    if hasattr(self,'_r_fine'):
      del self._r_fine
    self._r = 10*self._r

  def avg_y(self):
    '''
    Get the average of fiber coordinates along y axis

    Returns:
      float
    '''
    return np.mean(self.y())

  def x(self, mode='rough') -> np.ndarray:
    '''
    get x coordinates

    Parameter:
      mode (str): determines if 'rough' mesh is returned or 'fine' mesh is returned
    '''
    return self.r(mode)[:, 0]

  def y(self, mode='rough') -> np.ndarray:
    '''
    get y coordinates

    Parameter:
      mode (str): determines if 'rough' mesh is returned or 'fine' mesh is returned
    '''
    return self.r(mode)[:, 1]

  def z(self, mode='rough') -> np.ndarray:
    '''
    get z coordinates

    Parameter:
      mode (str): determines if 'rough' mesh is returned or 'fine' mesh is returned
    '''
    return self.r(mode)[:, 2]
