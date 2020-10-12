import numpy as np

def cartesian_basis_vector(z_axis):
  """
  returns a cartesian coordinate basis with z_axis as the vector along z_axis

  Parameters
  ----------
  z_axis : numpy.ndarray of size 3

  Returns
  -------
  a1 : numpy.ndarray of shape (1,3) along x-axis. a1 is normalized
  a2 : numpy.ndarray of shape (1,3) along y-axis. a2 is normalized
  a3 : numpy.ndarray of shape (1,3) along z-axis. a3 is normalized
  """
  assert (np.linalg.norm(z_axis) > 0), "a3 vector should not have magnitude of zero."
  assert (z_axis.size == 3), "a3 should have length of 3"
  
  a3 = np.reshape(z_axis, (1, 3))/np.linalg.norm(z_axis)

  a1 = np.zeros((1, 3))
  while(np.linalg.norm(a1) == 0):
      a1 = np.random.uniform(-1, 1, (1, 3))
      a1 -= a1.dot(a3.T)*a3
  a1 /= np.linalg.norm(a1)

  a2 = np.cross(a3, a1)
  a2 /= np.linalg.norm(a2)
  
  return a1, a2, a3

def direction_cosine(axis):
  """
  get direction cosine for an input vector
  input should be an array of length 3
  """
  assert (axis.size is 3), f"axis should have length of 3 but is: {axis.size}"
  
  n = np.linalg.norm(axis)
  assert (n > 0), 'axis vector cannot be zero'

  axis = axis / n

  cosine = np.zeros((3, 1))
  cosine[0] = axis.dot([1, 0, 0])
  cosine[1] = axis.dot([0, 1, 0])
  cosine[2] = axis.dot([0, 0, 1])
  return cosine

def get_quaternion(axis, alpha):
  """
  given and input axis and angle alpha, calculate quaternion corresponding to clockwise rotation
  around input axis when looking along positive direction of axis
  look here: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  """
  if alpha == 0:
    return np.array([1, 0, 0, 0])

  beta = direction_cosine(axis)
  quaternion = np.zeros((4, 1))
  quaternion[0] = np.cos(alpha/2)
  quaternion[1] = np.sin(alpha/2)*beta[0]
  quaternion[2] = np.sin(alpha/2)*beta[1]
  quaternion[3] = np.sin(alpha/2)*beta[2]
  return quaternion

def rotation_from_quaternion(q):
  """
  returns rotation matrix for an arbitrary rotation

  Returns
  -------
  q: calculate the rotation matrix from an input quaternion that describes a rotation

  Notes
  -----
  look here for more:
  https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  """
  R = np.zeros((3, 3))
  R[0, 0] = q[0]**2+q[1]**2-q[2]**2-q[3]**2
  R[0, 1] = 2*(q[1]*q[2]-q[0]*q[3])
  R[0, 2] = 2*(q[0]*q[2]+q[1]*q[3])
  R[1, 0] = 2*(q[1]*q[2]+q[0]*q[3])
  R[1, 1] = q[0]**2-q[1]**2+q[2]**2-q[3]**2
  R[1, 2] = 2*(q[2]*q[3]-q[0]*q[1])
  R[2, 0] = 2*(q[1]*q[3]-q[0]*q[2])
  R[2, 1] = 2*(q[0]*q[1]+q[2]*q[3])
  R[2, 2] = q[0]**2-q[1]**2-q[2]**2+q[3]**2
  return R

def HCP_coordinates(diameter=5, lattice_constant=1) -> np.ndarray:
  """
  Get coordinates of a hexagonal close-packed (HCP) lattice in a circular area.
  We do this by using a BFS algorithm to check all the neighbors.
  
  Parameters
  ----------
    diameter : diameter of the circular area
    lattice_constant : lattice constant for the HCP lattice

  Returns
  -------
    coordinates : numpy.ndarray with shape (N,2) where N is the number of points in the lattice
  """
  coordinates = np.zeros((0, 2))
  a = [np.array([1, 0]), np.array([np.cos(np.pi/3), np.sin(np.pi/3)])]
  a = [base*lattice_constant for base in a]

  s = [[0, 0]]
  visited = []
  while(len(s) > 0):
    node = s.pop()
    if (node in visited):
      continue

    c = node[0]*a[0]+node[1]*a[1]
    if (np.linalg.norm(c) < diameter):
      coordinates = np.vstack((coordinates, c))
      visited.append(node)

      nextNode = [node[0]+1, node[1]+0]
      if (not nextNode in visited):
        s.append(nextNode)

      nextNode = [node[0]+0, node[1]+1]
      if (not nextNode in visited):
        s.append(nextNode)

      nextNode = [node[0]-1, node[1]+1]
      if (not nextNode in visited):
        s.append(nextNode)

      nextNode = [node[0]-1, node[1]+0]
      if (not nextNode in visited):
        s.append(nextNode)

      nextNode = [node[0]+0, node[1]-1]
      if (not nextNode in visited):
        s.append(nextNode)

      nextNode = [node[0]+1, node[1]-1]
      if (not nextNode in visited):
        s.append(nextNode)

  return coordinates
