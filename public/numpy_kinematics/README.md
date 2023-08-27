# Rosetta-Style Kinematics in Numpy Syntax

KEYWORDS: CORE_CONCEPTS DOCKING

As toolchains evolve, it is useful to compare a set of unit tests against their libraries to make sure we can reproduce Rosetta's logic. These unit tests were created by calling `core/kinematics/*` functions and objects, but operate in Numpy syntax.

This uses flax's numpy-like syntax for the purpose of unit testing. The code should still work if you just use `import numpy as np` and leave out some jax/flax bits.

```py
import jax
import jax.numpy as jnp
import jax.numpy as np
import flax.linen as nn

class BatchDotProduct(nn.Module):
    @nn.compact
    def __call__(self, x: jnp.ndarray, y: jnp.ndarray) -> jnp.ndarray:
        return np.sum(x*y, axis=-1, keepdims=True)
        pass

    @classmethod
    def debug(cls):
        a = np.array([[1,2,3], [3,4,5]])
        b = np.array([[1,2,3], [1,2,3]])
        o = np.array([ [14.],   [26.], ])
        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), a, b )
        output = layer.apply( variables, a, b )
        print( output )
        assert jnp.allclose(
            output,
            o
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class Angle(nn.Module):
    @nn.compact
    def __call__(self,
                 p0: jnp.ndarray,
                 p1: jnp.ndarray,
                 p2: jnp.ndarray,
                 ) -> jnp.ndarray:
        a = p1 - p0
        b = p1 - p2
        term1 = np.linalg.norm(
            np.cross(a,b),
            axis=-1,
            keepdims=True
        )
        term2 = BatchDotProduct()(a,b)
        return np.abs( np.arctan2(term2, term1) )

    @classmethod
    def debug(cls):
        p0 = jnp.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p1 = jnp.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p2 = jnp.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), p0, p1, p2 )
        output = layer.apply( variables, p0, p1, p2 )
        print( output )

        assert jnp.allclose(
            output, 
            jnp.array([[0.44907826],
                       [0.5157146 ],
                       [0.4016324 ],
                       [0.36969098]], dtype='float32')
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class Dihedral(nn.Module):
    @nn.compact
    def __call__(self,
                 p0: jnp.ndarray,
                 p1: jnp.ndarray,
                 p2: jnp.ndarray,
                 p3: jnp.ndarray
                 ) -> jnp.ndarray:
        # https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        # normalize b1 so that it does not influence magnitude of vector
        # rejections that come next
        b1 /= np.linalg.norm(b1, axis=-1, keepdims=True)

        # vector rejections
        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - BatchDotProduct()(b0, b1)*b1
        w = b2 - BatchDotProduct()(b2, b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = BatchDotProduct()(v, w)
        y = BatchDotProduct()(np.cross(b1, v), w)
        return np.arctan2(y, x)

    @classmethod
    def debug(cls):
        p0 = jnp.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])
        p1 = jnp.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p2 = jnp.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p3 = jnp.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), p0, p1, p2, p3 )
        output = layer.apply( variables, p0, p1, p2, p3 )
        print( output )

        assert jnp.allclose(
            output, 
            jnp.array([[ 1.8286608],
                       [-0.503368 ],
                       [ 1.4385247],
                       [ 1.8120883]],
                      dtype='float32')
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class ComputeInputStub(nn.Module):
    @nn.compact
    def __call__(self, ggp: jnp.ndarray, gp: jnp.ndarray, p: jnp.ndarray ):
        """
        Computes the stub used to create an atom from its ancestors. Mirrors Rosetta methodology to a degree.

        Parameters
        ----------
        ggp : jnp.ndarray
            XYZ coords of shape (num_nodes, 3) of the great-grandparent atom.
        gp : jnp.ndarray
            XYZ coords of shape (num_nodes, 3) of the grandparent atom.
        p : jnp.ndarray
            XYZ coords of shape (num_nodes, 3) of the parent atom.

        Returns
        -------
        M : np.ndarray
            Rotation Matrix of shape (num_nodes, 3, 3).
        v : np.ndarray
            Rotation Matrix of shape (num_nodes, 3).
        """
        a = p
        b = gp
        c = ggp
        center = p

        e1 = a-b
        e1 = e1 / np.linalg.norm( e1, axis=-1, keepdims=True )

        e3 = np.cross( e1, c-b )
        e3 = e3 / np.linalg.norm( e3, axis=-1, keepdims=True )

        e2 = np.cross( e3, e1 )

        e1 = np.expand_dims(e1, -1)
        e2 = np.expand_dims(e2, -1)
        e3 = np.expand_dims(e3, -1)

        return np.concatenate( [e1,e2,e3], axis=-1 ), p


    #@staticmethod
    @classmethod
    def debug(cls):
        # first, santity check that we know how to calc distances
        a = jnp.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2]
        ])
        b = jnp.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2]
        ])
        diff = a-b
        dist = np.linalg.norm( diff, axis=-1, keepdims=True )
        assert dist[0][0] > 0.99 # 1
        assert dist[0][0] < 1.01 # 1
        assert dist[1][0] > 4.99 # 5
        assert dist[1][0] < 5.01 # 5
        assert dist[2][0] > 6.23 # 6.23298
        assert dist[2][0] < 6.24 # 6.23298
        print( dist )

        ggp = a
        gp = b
        p = jnp.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2]
        ])


        layer = ComputeInputStub()
        variables = layer.init( jax.random.PRNGKey(0), ggp, gp, p )
        output = layer.apply( variables, ggp, gp, p )

        # v
        assert jnp.allclose(
            output[1], 
            jnp.array([[20.  , 12.  , 12.  ],
                       [22.  ,  2.  , 42.  ],
                       [ 2.21,  2.2 , -2.2 ]], dtype='float32')
        )

        # M
        assert jnp.allclose(
            output[0], 
            jnp.array([[[ 0.78933704,  0.3803963 , -0.48191875],
                        [ 0.43413538,  0.20921798,  0.8762159 ],
                        [ 0.43413538, -0.9008476 ,  0.        ]],

                       [[ 0.46524212,  0.26373896, -0.8449802 ],
                        [ 0.04652421,  0.945977  ,  0.32087854],
                        [ 0.88396   , -0.18859825,  0.42783806]],

                       [[ 0.04913414,  0.06628565, -0.99659026],
                        [-0.44667345, -0.89099705, -0.0812844 ],
                        [-0.8933469 ,  0.44914424, -0.01417033]]], dtype='float32')
        )
        print( "passed {} tests".format( cls.__name__ ) )
        pass


class XRotMatRadians(nn.Module):
    @nn.compact
    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        #x = np.expand_dims(x, -2)
        angle_sin = np.sin( x )
        angle_cos = np.cos( x )
    
        ones = np.ones_like( x )
        zeros = np.zeros_like( x )

        angle_sin_neg = -1 * angle_sin

        # doing columns instead of rows, even though the code reads like rows
        col_x = np.concatenate( [ones,zeros,zeros], axis=-1 )
        #print( col_x.shape )
        col_y = np.concatenate( [zeros,angle_cos,angle_sin], axis=-1 )
        col_z = np.concatenate( [zeros,angle_sin_neg,angle_cos], axis=-1 )

        #col_x = np.expand_dims(col_x, -1)
        #col_y = np.expand_dims(col_y, -1)
        #col_z = np.expand_dims(col_z, -1)
    
        return np.stack( [col_x,col_y,col_z], axis=-1 )


    @classmethod
    def debug(cls):
        a = jnp.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), a )
        output = layer.apply( variables, a )
        assert jnp.allclose(
            output, 
            jnp.array([[[ 1.       , 0.       , 0.        ],
                        [ 0.       , 0.5403023,-0.84147096],
                        [ 0.       , 0.84147096,  0.5403023 ]],
                       
                       [[ 1.       , 0.       , 0.        ],
                        [ 0.       ,-0.41614684, -0.9092974 ],
                        [ 0.       , 0.9092974,-0.41614684]],

                       [[ 1.       , 0.       , 0.        ],
                        [ 0.       ,-0.12884454, -0.9916648 ],
                        [ 0.       , 0.9916648,-0.12884454]],
                       
                       [[ 1.       , 0.       , 0.        ],
                        [ 0.       , 1.       ,-0.        ],
                        [ 0.       , 0.       , 1.        ]],
                       
                       [[ 1.       , 0.       , 0.        ],
                        [ 0.       ,-0.61248755,  0.7904802 ],
                        [ 0.       ,-0.7904802,-0.61248755]]],
                      dtype='float32')

        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class ZRotMatRadians(nn.Module):
    @nn.compact
    def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
        #x = np.expand_dims(x, -2)
        angle_sin = np.sin( x )
        angle_cos = np.cos( x )
    
        ones = np.ones_like( x )
        zeros = np.zeros_like( x )

        angle_sin_neg = -1 * angle_sin

        # doing columns instead of rows, even though the code reads like rows
        col_x = np.concatenate( [angle_cos,angle_sin,zeros], axis=-1 )
        col_y = np.concatenate( [angle_sin_neg,angle_cos,zeros], axis=-1 )
        col_z = np.concatenate( [zeros,zeros,ones], axis=-1 )

        #col_x = np.expand_dims(col_x, -1)
        #col_y = np.expand_dims(col_y, -1)
        #col_z = np.expand_dims(col_z, -1)
    
        return np.stack( [col_x,col_y,col_z], axis=-1 )


    @classmethod
    def debug(cls):
        a = jnp.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), a )
        output = layer.apply( variables, a )
        assert jnp.allclose(
            output, 
            jnp.array([[[ 0.5403023 , -0.84147096,  0.        ],
                        [ 0.84147096,  0.5403023 ,  0.        ],
                        [ 0.        ,  0.        ,  1.        ]],

                       [[-0.41614684, -0.9092974 ,  0.        ],
                        [ 0.9092974 , -0.41614684,  0.        ],
                        [ 0.        ,  0.        ,  1.        ]],

                       [[-0.12884454, -0.9916648 ,  0.        ],
                        [ 0.9916648 , -0.12884454,  0.        ],
                        [ 0.        ,  0.        ,  1.        ]],

                       [[ 1.        , -0.        ,  0.        ],
                        [ 0.        ,  1.        ,  0.        ],
                        [ 0.        ,  0.        ,  1.        ]],

                       [[-0.61248755,  0.7904802 ,  0.        ],
                        [-0.7904802 , -0.61248755,  0.        ],
                        [ 0.        ,  0.        ,  1.        ]]],
                      dtype='float32')

        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class TranslateFromStubB(nn.Module):
    @nn.compact
    def __call__(self, StubB: jnp.ndarray, ic_in: jnp.ndarray) -> jnp.ndarray:
        col_x = StubB[:,:,0]
        d = ic_in[:,2]
        d = np.expand_dims(d, -1)
        return d * col_x

    @classmethod
    def debug(cls):
        # Currently tested in GetAtomXYZ.debug
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class GetAtomXYZ(nn.Module):
    @nn.compact
    def __call__(self,
                 p: jnp.ndarray,
                 gp: jnp.ndarray,
                 ggp: jnp.ndarray,
                 ic_in: jnp.ndarray
                 ) -> jnp.ndarray:
        parent_M, parent_v = ComputeInputStub()(ggp,gp,p)

        phi = ic_in[:,0]
        phi = np.expand_dims(phi, -1)
        x_phi = XRotMatRadians()(phi)

        theta = ic_in[:,1]
        theta = np.expand_dims(theta, -1)
        z_theta = ZRotMatRadians()(theta)

        StubA = np.matmul( parent_M, x_phi )
        StubB = np.matmul( StubA, z_theta )

        dv = TranslateFromStubB()( StubB, ic_in )

        return parent_v + dv

    @classmethod
    def debug(cls):
        ggp = jnp.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        gp = jnp.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p = jnp.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])
        ic = jnp.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), p, gp, ggp, ic )
        output = layer.apply( variables, p, gp, ggp, ic )
        #print( output )

        assert jnp.allclose(
            output, 
            jnp.array([[20.07297  , 12.033433 , 12.059646 ],
                       [24.958748 ,  1.6791693, 51.54688  ],
                       [ 2.1709127,  2.0679247, -2.0549886],
                       [ 2.2489984,  2.3320887, -2.445023 ]],
                      dtype='float32')
        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class StubFromCoords(nn.Module):
    @nn.compact
    def __call__(self,
                 N: jnp.ndarray,
                 CA: jnp.ndarray,
                 C: jnp.ndarray
                 ) -> jnp.ndarray:

        #matching vocab from Stub.cc
        a = N
        b = CA
        c = C
        center = CA

        e1 = a-b
        e1 = e1 / np.linalg.norm( e1, axis=-1, keepdims=True )

        e3 = np.cross( e1, c-b )
        e3 = e3 / np.linalg.norm( e3, axis=-1, keepdims=True )

        e2 = np.cross( e3, e1 )

        e1 = np.expand_dims(e1, -1)
        e2 = np.expand_dims(e2, -1)
        e3 = np.expand_dims(e3, -1)

        #foo = np.concatenate( [e1,e2,e3], axis=-1 )
        #center = np.expand_dims( center, -1 )
        #foo = np.concatenate( [foo,center], axis=-1 )
        #print( foo.shape )

        return np.concatenate( [e1,e2,e3], axis=-1 ), center

        pass

    @classmethod
    def debug(cls):
        N = jnp.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        CA = jnp.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        C = jnp.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        o0 = jnp.array([[[ 0.        ,  0.87621593,  0.48191875],
                         [ 0.        ,  0.48191875, -0.87621593],
                         [-1.        ,  0.        ,  0.        ]],

                        [[ 0.        ,  0.5347976 ,  0.8449802 ],
                         [ 0.8       ,  0.5069881 , -0.32087854],
                         [-0.6       ,  0.67598414, -0.42783806]],

                        [[ 0.08021849,  0.01931177,  0.99659014],
                         [-0.99470925, -0.06281925,  0.08128439],
                         [ 0.06417479, -0.99783796,  0.01417033]],

                        [[ 0.0800128 ,  0.02168863,  0.9965579 ],
                         [-0.9921587 , -0.09456854,  0.08171775],
                         [ 0.09601536, -0.99528205,  0.01395187]]], dtype='float32')

        o1 = jnp.array([[ 0. ,  1. ,  1. ],
                        [ 2. ,  0. ,  4. ],
                        [ 2.1,  3.2, -0.2],
                        [ 2.1,  3.2, -0.3]], dtype='float32')

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), N, CA, C )
        output = layer.apply( variables, N, CA, C )
        #print( output )

        assert jnp.allclose( output[0], o0 )
        assert jnp.allclose( output[1], o1 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class HTValsFromCoords(nn.Module):
    @nn.compact
    def __call__(self,
                 Ni: jnp.ndarray,
                 CAi: jnp.ndarray,
                 Ci: jnp.ndarray,
                 Nf: jnp.ndarray,
                 CAf: jnp.ndarray,
                 Cf: jnp.ndarray
                 ) -> jnp.ndarray:
        Mi, Vi = StubFromCoords()(Ni, CAi, Ci)
        Mf, Vf = StubFromCoords()(Nf, CAf, Cf)

        Mi_t = np.transpose( Mi, axes=[0,2,1] )

        dV = Vf - Vi
        dV = np.expand_dims(dV, -1)

        T = np.matmul( Mi_t, dV )
        R = np.matmul( Mi_t, Mf )

        CAdist = np.linalg.norm(CAi-CAf, axis=-1, keepdims=True)
        #print( T.shape, R.shape, CAdist.shape )
        # (4, 3, 1) (4, 3, 3) (4, 1)
        #exit( 0 )

        Tflat = np.reshape( T, newshape=[-1, 3] )
        Rflat = np.reshape( R, newshape=[-1, 9] )

        return np.concatenate( [Tflat,Rflat,CAdist], axis=-1 )

    @classmethod
    def debug(cls):
        p0 = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p1 = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p2 = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])
        p3 = np.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ]) 
        o = jnp.array([[-1.1000000e+01,  2.2825424e+01,  0.0000000e+00,  4.3413538e-01,
                        4.4217414e-01, -7.8486210e-01, -9.0084767e-01,  2.1309203e-01,
                        -3.7823996e-01,  2.9802322e-08,  8.7124860e-01,  4.9084231e-01,
                        2.5337719e+01],
                       [-2.1200001e+01,  3.7397324e+01, -1.9073486e-06,  4.9315664e-01,
                        -2.9688621e-01,  8.1771332e-01, -8.6994052e-01, -1.6830048e-01,
                        4.6354976e-01, -2.9802322e-08, -9.3996465e-01, -3.4127179e-01,
                        4.2988369e+01],
                       [ 8.7518364e-01,  2.0606194e+00,  3.7252903e-09, -3.9092135e-01,
                         7.8699952e-01, -4.7729683e-01, -9.2042404e-01, -3.3425343e-01,
                         2.0271687e-01, -1.8626451e-09,  5.1856184e-01,  8.5504007e-01,
                         2.2387719e+00],
                       [ 8.0892938e-01,  2.0875185e+00,  3.7252903e-09, -3.6132732e-01,
                         8.1903273e-01, -4.4567707e-01, -9.3243903e-01, -3.1738147e-01,
                         1.7270330e-01, -9.3132257e-10,  4.7796917e-01,  8.7837678e-01,
                         2.2387719e+00]], dtype='float32')

        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), p0, p1, p2, p1, p2, p3 )
        output = layer.apply( variables, p0, p1, p2, p1, p2, p3 )
        #print( output )

        assert jnp.allclose( output, o, atol=1e-04 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class BuildHTs(nn.Module):
    @nn.compact
    def __call__(self,
                 A: jnp.ndarray,
                 N: jnp.ndarray,
                 CA: jnp.ndarray,
                 C: jnp.ndarray
                 ) -> jnp.ndarray:
        index_sources = A[:,0]
        index_targets = A[:,1]
        # from jax.lax.gather:
        # The semantics of gather are complicated, and its API might change in the future. For most use cases, you should prefer Numpy-style indexing (e.g., x[:, (1,4,7), ...]), rather than using gather directly.

        Nsource = N[index_sources,:]
        CAsource = CA[index_sources,:]
        Csource = C[index_sources,:]

        Ntarget = N[index_targets,:]
        CAtarget = CA[index_targets,:]
        Ctarget = C[index_targets,:]
        
        return HTValsFromCoords()(Nsource,CAsource,Csource,Ntarget,CAtarget,Ctarget)

    @classmethod
    def debug(cls):
        N = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        CA = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        C = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        A = np.array([
            [0, 1],
            [1, 0],
            [1, 2],
            [3, 1]
        ])

        Ni = np.array([
            [0., 1., 0.], #0
            [2., 4., 1.,], #1
            [2., 4., 1.,], #1
            [2.6, -3.0, 0.3], #3
        ])

        Nf = np.array([
            [2., 4., 1.,], #1
            [0., 1., 0.], #0
            [2.6, -3.0, 0.2], #2
            [2., 4., 1.,], #1
        ])

        CAi = np.array([
            [0., 1., 1.], #0
            [2., 0., 4.,], #1
            [2., 0., 4.,], #1
            [2.1, 3.2, -0.3] #3
        ])

        CAf = np.array([
            [2., 0., 4.,], #1
            [0., 1., 1.], #0
            [2.1, 3.2, -0.2], #2
            [2., 0., 4.,], #1
        ])

        Ci = np.array([
            [20., 12., 12.], #0
            [22., 2., 42.,], #1
            [22., 2., 42.,], #1
            [2.21, 2.2, -2.3] #3
        ])

        Cf = np.array([
            [22., 2., 42.,], #1
            [20., 12., 12.], #0
            [2.21, 2.2, -2.2], #2
            [22., 2., 42.,], #1
        ])


        layer = cls()
        variables = layer.init( jax.random.PRNGKey(0), A, N, CA, C )
        output = layer.apply( variables, A, N, CA, C )

        layer2 = HTValsFromCoords()
        var2 = layer2.init( jax.random.PRNGKey(0), Ni, CAi, Ci, Nf, CAf, Cf )
        o2 = layer2.apply( var2, Ni, CAi, Ci, Nf, CAf, Cf )

        assert jnp.allclose( output, o2 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

if __name__ == '__main__':
    BatchDotProduct.debug()
    Angle.debug()
    Dihedral.debug()
    ComputeInputStub.debug()
    XRotMatRadians.debug()
    ZRotMatRadians.debug()
    TranslateFromStubB.debug()
    GetAtomXYZ.debug()
    StubFromCoords.debug()
    HTValsFromCoords.debug()
    BuildHTs.debug()
```