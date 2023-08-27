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

## PyTorch

```py
import numpy as np
import torch
from torch import nn

class BatchDotProduct(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x, y):
        return torch.sum(x*y, dim=-1, keepdim=True)

    @classmethod
    def debug(cls):
        a = np.array([[1,2,3], [3,4,5]])
        b = np.array([[1,2,3], [1,2,3]])
        o = np.array([ [14.],   [26.], ])

        a = torch.from_numpy(a)
        b = torch.from_numpy(b)

        layer = cls()
        output = layer( a, b ).numpy()

        np.testing.assert_array_almost_equal(
            output,
            o
        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class Angle(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, p0, p1, p2):
        a = p1 - p0
        b = p1 - p2
        term1 = torch.linalg.norm(
            torch.linalg.cross(a,b),
            dim=-1,
            keepdims=True
        )
        term2 = BatchDotProduct()(a,b)
        return torch.abs( torch.arctan2(term2, term1) )

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

        p0 = torch.from_numpy(p0)
        p1 = torch.from_numpy(p1)
        p2 = torch.from_numpy(p2)

        layer = cls()
        output = layer( p0, p1, p2 ).numpy()
        print( output )

        np.testing.assert_array_almost_equal(
            output, 
            np.array([[0.44907826],
                       [0.5157146 ],
                       [0.4016324 ],
                       [0.36969098]], dtype='float32')
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class Dihedral(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, p0, p1, p2, p3 ):
        # https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        # normalize b1 so that it does not influence magnitude of vector
        # rejections that come next
        b1 /= torch.linalg.norm(b1, dim=-1, keepdims=True)

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
        y = BatchDotProduct()(torch.linalg.cross(b1, v), w)
        return torch.arctan2(y, x)

    @classmethod
    def debug(cls):
        p0 = np.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])
        p1 = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p2 = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p3 = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        p0 = torch.from_numpy(p0)
        p1 = torch.from_numpy(p1)
        p2 = torch.from_numpy(p2)
        p3 = torch.from_numpy(p3)

        layer = cls()
        output = layer( p0, p1, p2, p3 ).numpy()
        print( output )

        np.testing.assert_array_almost_equal(
            output, 
            np.array([[ 1.8286608],
                       [-0.503368 ],
                       [ 1.4385247],
                       [ 1.8120883]],
                      dtype='float32')
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass


class ComputeInputStub(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, ggp, gp, p ):
        """
        Computes the stub used to create an atom from its ancestors. Mirrors Rosetta methodology to a degree.

        Parameters
        ----------
        ggp 
            XYZ coords of shape (num_nodes, 3) of the great-grandparent atom.
        gp 
            XYZ coords of shape (num_nodes, 3) of the grandparent atom.
        p 
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
        e1 = e1 / torch.linalg.norm( e1, dim=-1, keepdims=True )

        e3 = torch.linalg.cross( e1, c-b )
        e3 = e3 / torch.linalg.norm( e3, dim=-1, keepdims=True )

        e2 = torch.linalg.cross( e3, e1 )

        e1 = torch.unsqueeze(e1, -1)
        e2 = torch.unsqueeze(e2, -1)
        e3 = torch.unsqueeze(e3, -1)

        return torch.concat( [e1,e2,e3], dim=-1 ), p


    #@staticmethod
    @classmethod
    def debug(cls):
        # first, santity check that we know how to calc distances
        a = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2]
        ])
        b = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2]
        ])
        a = torch.from_numpy(a)
        b = torch.from_numpy(b)

        diff = a-b
        dist = torch.linalg.norm( diff, dim=-1, keepdims=True )
        assert dist[0][0] > 0.99 # 1
        assert dist[0][0] < 1.01 # 1
        assert dist[1][0] > 4.99 # 5
        assert dist[1][0] < 5.01 # 5
        assert dist[2][0] > 6.23 # 6.23298
        assert dist[2][0] < 6.24 # 6.23298
        print( dist )

        ggp = a
        gp = b
        p = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2]
        ])

        p = torch.from_numpy(p)

        layer = ComputeInputStub()
        output = layer( ggp, gp, p )

        # v
        np.testing.assert_array_almost_equal(
            output[1], 
            np.array([[20.  , 12.  , 12.  ],
                       [22.  ,  2.  , 42.  ],
                       [ 2.21,  2.2 , -2.2 ]], dtype='float32')
        )

        # M
        np.testing.assert_array_almost_equal(
            output[0], 
            np.array([[[ 0.78933704,  0.3803963 , -0.48191875],
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

    def __init__(self):
        super().__init__()

    def forward(self, x):
        #x = torch.unsqueeze(x, -2)
        angle_sin = torch.sin( x )
        angle_cos = torch.cos( x )
    
        ones = torch.ones_like( x )
        zeros = torch.zeros_like( x )

        angle_sin_neg = -1 * angle_sin

        # doing columns instead of rows, even though the code reads like rows
        col_x = torch.concat( [ones,zeros,zeros], dim=-1 )
        #print( col_x.shape )
        col_y = torch.concat( [zeros,angle_cos,angle_sin], dim=-1 )
        col_z = torch.concat( [zeros,angle_sin_neg,angle_cos], dim=-1 )

        #col_x = torch.unsqueeze(col_x, -1)
        #col_y = torch.unsqueeze(col_y, -1)
        #col_z = torch.unsqueeze(col_z, -1)
    
        return torch.stack( [col_x,col_y,col_z], dim=-1 )


    @classmethod
    def debug(cls):
        a = np.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])
        a = torch.from_numpy(a)

        layer = cls()
        output = layer( a ).numpy()
        np.testing.assert_array_almost_equal(
            output, 
            np.array([[[ 1.       , 0.       , 0.        ],
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

    def __init__(self):
        super().__init__()

    def forward(self, x):
        #x = torch.unsqueeze(x, -2)
        angle_sin = torch.sin( x )
        angle_cos = torch.cos( x )
    
        ones = torch.ones_like( x )
        zeros = torch.zeros_like( x )

        angle_sin_neg = -1 * angle_sin

        # doing columns instead of rows, even though the code reads like rows
        col_x = torch.concat( [angle_cos,angle_sin,zeros], dim=-1 )
        col_y = torch.concat( [angle_sin_neg,angle_cos,zeros], dim=-1 )
        col_z = torch.concat( [zeros,zeros,ones], dim=-1 )

        #col_x = torch.unsqueeze(col_x, -1)
        #col_y = torch.unsqueeze(col_y, -1)
        #col_z = torch.unsqueeze(col_z, -1)
    
        return torch.stack( [col_x,col_y,col_z], dim=-1 )


    @classmethod
    def debug(cls):
        a = np.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])
        a = torch.from_numpy(a)

        layer = cls()
        output = layer( a ).numpy()
        np.testing.assert_array_almost_equal(
            output, 
            np.array([[[ 0.5403023 , -0.84147096,  0.        ],
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

    def __init__(self):
        super().__init__()

    def forward(self, StubB, ic_in):
        col_x = StubB[:,:,0]
        d = ic_in[:,2]
        d = torch.unsqueeze(d, -1)
        return d * col_x

    @classmethod
    def debug(cls):
        # Currently tested in GetAtomXYZ.debug
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class GetAtomXYZ(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, p, gp, ggp, ic_in ):
        parent_M, parent_v = ComputeInputStub()(ggp,gp,p)

        phi = ic_in[:,0]
        phi = torch.unsqueeze(phi, -1)
        x_phi = XRotMatRadians()(phi)

        theta = ic_in[:,1]
        theta = torch.unsqueeze(theta, -1)
        z_theta = ZRotMatRadians()(theta)

        StubA = torch.matmul( parent_M, x_phi )
        StubB = torch.matmul( StubA, z_theta )

        dv = TranslateFromStubB()( StubB, ic_in )

        return parent_v + dv

    @classmethod
    def debug(cls):
        ggp = torch.tensor([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        gp = torch.tensor([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p = torch.tensor([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])
        ic = torch.tensor([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])

        layer = cls()
        output = layer( p, gp, ggp, ic ).numpy()

        np.testing.assert_array_almost_equal(
            output, 
            np.array([[20.07297  , 12.033433 , 12.059646 ],
                       [24.958748 ,  1.6791693, 51.54688  ],
                       [ 2.1709127,  2.0679247, -2.0549886],
                       [ 2.2489984,  2.3320887, -2.445023 ]],
                      dtype='float32')
        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass


class StubFromCoords(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self,
                 N,
                 CA,
                 C
                 ):

        #matching vocab from Stub.cc
        a = N
        b = CA
        c = C
        center = CA

        e1 = a-b
        e1 = e1 / torch.linalg.norm( e1, dim=-1, keepdims=True )

        e3 = torch.linalg.cross( e1, c-b )
        e3 = e3 / torch.linalg.norm( e3, dim=-1, keepdims=True )

        e2 = torch.linalg.cross( e3, e1 )

        e1 = torch.unsqueeze(e1, -1)
        e2 = torch.unsqueeze(e2, -1)
        e3 = torch.unsqueeze(e3, -1)

        #foo = torch.concat( [e1,e2,e3], dim=-1 )
        #center = torch.unsqueeze( center, -1 )
        #foo = torch.concat( [foo,center], dim=-1 )
        #print( foo.shape )

        return torch.concat( [e1,e2,e3], dim=-1 ), center

        pass

    @classmethod
    def debug(cls):
        N = torch.tensor([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        CA = torch.tensor([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        C = torch.tensor([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        o0 = np.array([[[ 0.        ,  0.87621593,  0.48191875],
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

        o1 = np.array([[ 0. ,  1. ,  1. ],
                        [ 2. ,  0. ,  4. ],
                        [ 2.1,  3.2, -0.2],
                        [ 2.1,  3.2, -0.3]], dtype='float32')

        layer = cls()
        output = [l.numpy() for l in layer( N, CA, C )]

        np.testing.assert_array_almost_equal( output[0], o0 )
        np.testing.assert_array_almost_equal( output[1], o1 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass


class HTValsFromCoords(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self,
                Ni, CAi, Ci,
                Nf, CAf, Cf ):
        Mi, Vi = StubFromCoords()(Ni, CAi, Ci)
        Mf, Vf = StubFromCoords()(Nf, CAf, Cf)

        Mi_t = torch.transpose( Mi, dim0=1, dim1=2 )

        dV = Vf - Vi
        dV = torch.unsqueeze(dV, -1)

        T = torch.matmul( Mi_t, dV )
        R = torch.matmul( Mi_t, Mf )

        CAdist = torch.linalg.norm(CAi-CAf, dim=-1, keepdims=True)
        #print( T.shape, R.shape, CAdist.shape )
        # (4, 3, 1) (4, 3, 3) (4, 1)
        #exit( 0 )

        Tflat = torch.reshape( T, shape=[-1, 3] )
        Rflat = torch.reshape( R, shape=[-1, 9] )

        return torch.concat( [Tflat,Rflat,CAdist], dim=-1 )

    @classmethod
    def debug(cls):
        p0 = torch.tensor([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p1 = torch.tensor([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p2 = torch.tensor([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])
        p3 = torch.tensor([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ]) 
        o = np.array([[-1.1000000e+01,  2.2825424e+01,  0.0000000e+00,  4.3413538e-01,
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
        output = layer( p0, p1, p2, p1, p2, p3 ).numpy()
        #print( output )

        np.testing.assert_array_almost_equal( output, o )

        print( "passed {} tests".format( cls.__name__ ) )
        pass


class BuildHTs(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, A, N, CA, C ):
        index_sources = A[:,0]
        index_targets = A[:,1]

        Nsource = N[index_sources,:]
        CAsource = CA[index_sources,:]
        Csource = C[index_sources,:]

        Ntarget = N[index_targets,:]
        CAtarget = CA[index_targets,:]
        Ctarget = C[index_targets,:]
        
        return HTValsFromCoords()(Nsource,CAsource,Csource,Ntarget,CAtarget,Ctarget)

    @classmethod
    def debug(cls):
        N = torch.tensor([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        CA = torch.tensor([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        C = torch.tensor([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        A = torch.tensor([
            [0, 1],
            [1, 0],
            [1, 2],
            [3, 1]
        ])

        Ni = torch.tensor([
            [0., 1., 0.], #0
            [2., 4., 1.,], #1
            [2., 4., 1.,], #1
            [2.6, -3.0, 0.3], #3
        ])

        Nf = torch.tensor([
            [2., 4., 1.,], #1
            [0., 1., 0.], #0
            [2.6, -3.0, 0.2], #2
            [2., 4., 1.,], #1
        ])

        CAi = torch.tensor([
            [0., 1., 1.], #0
            [2., 0., 4.,], #1
            [2., 0., 4.,], #1
            [2.1, 3.2, -0.3] #3
        ])

        CAf = torch.tensor([
            [2., 0., 4.,], #1
            [0., 1., 1.], #0
            [2.1, 3.2, -0.2], #2
            [2., 0., 4.,], #1
        ])

        Ci = torch.tensor([
            [20., 12., 12.], #0
            [22., 2., 42.,], #1
            [22., 2., 42.,], #1
            [2.21, 2.2, -2.3] #3
        ])

        Cf = torch.tensor([
            [22., 2., 42.,], #1
            [20., 12., 12.], #0
            [2.21, 2.2, -2.2], #2
            [22., 2., 42.,], #1
        ])


        layer = cls()
        output = layer( A, N, CA, C ).numpy()

        layer2 = HTValsFromCoords()
        o2 = layer2( Ni, CAi, Ci, Nf, CAf, Cf ).numpy()

        np.testing.assert_array_almost_equal( output, o2 )

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

## Tensorflow

```py
import numpy as np

import keras
from keras.layers import *
from keras.models import Model

from keras import backend
from keras import constraints
from keras import initializers
from keras import regularizers
from keras.engine.base_layer import Layer
from keras.engine.input_spec import InputSpec
from keras.utils import tf_utils

import tensorflow as tf

from tensorflow.python.framework import dtypes
from tensorflow.python.ops import math_ops
from tensorflow.python.util.tf_export import keras_export

class Angle(Layer):
  # PROPER bond angles
  # GetAtomXYZ uses IMPROPER bond angles
  
  def __init__(self,
               **kwargs):
    super(Angle, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    p1 = inputs[0] #(batch, 3)
    p2 = inputs[1] #(batch, 3)
    p3 = inputs[2] #(batch, 3)

    # https://www.mathworks.com/matlabcentral/answers/16243-angle-between-two-vectors-in-3d
    # https://stackoverflow.com/questions/21483999/using-atan2-to-find-angle-between-two-vectors
    a = p2 - p1
    b = p2 - p3
    
    term1 = tf.linalg.norm(
      tf.linalg.cross(a,b),
      axis=-1,
      keepdims=True
    )
    term2 = Dot(axes=-1)([a,b])
    return tf.math.abs( tf.math.atan2( x=term1, y=term2 ) )

  def get_config(self):
    config = {}
    base_config = super(Angle, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(1)

  @classmethod
  def debug(cls):
    #pylint: disable = E, W, R, C
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

    goal = np.array([[0.44907826],
             [0.5157146 ],
             [0.4016324 ],
             [0.36969098]], dtype='float32')
    
    layer = cls()
    output = layer([ p0, p1, p2 ]).numpy()

    np.testing.assert_array_almost_equal( output, goal )
    
    print( "passed {} tests".format( cls.__name__ ) )
    pass

class Dihedral(Layer):
  def __init__(self,
               apply_atan2 = True,
               **kwargs):
    super(Dihedral, self).__init__(**kwargs)
    self.apply_atan2 = apply_atan2
    
  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    p1 = inputs[0] #(batch, 3)
    p2 = inputs[1] #(batch, 3)
    p3 = inputs[2] #(batch, 3)
    p4 = inputs[3] #(batch, 3)

    a = p2 - p1
    b = p3 - p2
    c = p4 - p3

    a, _ = tf.linalg.normalize( a, axis=-1 )
    b, _ = tf.linalg.normalize( b, axis=-1 )
    c, _ = tf.linalg.normalize( c, axis=-1 )

    #print( a.shape )
    x_1 = Dot(axes=-1)([a,c])
    x_1 = tf.math.negative( x_1 )
    x_2 = Dot(axes=-1)([a,b])
    x_3 = Dot(axes=-1)([b,c])
    x = x_1 + ( x_2 * x_3 )

    y = Dot(axes=-1)([a,tf.linalg.cross( b, c )])
    if self.apply_atan2:
      return tf.math.atan2( y=y, x=x )
    else:
      return tf.concat( [x,y], axis=-1 )

  def get_config(self):
    config = {
      "apply_atan2" : self.apply_atan2
    }
    base_config = super(Dihedral, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(1)

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
        p0 = np.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])
        p1 = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        p2 = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p3 = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])

        layer = cls()
        output = layer([ p0, p1, p2, p3 ]).numpy()
        print( output )

        np.testing.assert_array_almost_equal(
            output, 
            np.array([[ 1.8286608],
                       [-0.503368 ],
                       [ 1.4385247],
                       [ 1.8120883]],
                      dtype='float32')
        )
        
        print( "passed {} tests".format( cls.__name__ ) )
        pass

class ComputeInputStub(Layer):
  def __init__(self,
               **kwargs):
    super(ComputeInputStub, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      p   = inputs[2] #position AND first axis
      gp  = inputs[1]
      ggp = inputs[0]

      #matching vocab from Stub.cc
      a = p
      b = gp
      c = ggp
      center = p

      e1 = a-b;
      e1, _ = tf.linalg.normalize( e1, axis=-1 )

      e3 = tf.linalg.cross( e1, c-b )
      e3, _ = tf.linalg.normalize( e3, axis=-1 )

      e2 = tf.linalg.cross( e3, e1 )

      e1 = tf.expand_dims(e1, -1)
      e2 = tf.expand_dims(e2, -1)
      e3 = tf.expand_dims(e3, -1)

      return tf.concat( [e1,e2,e3], axis=-1 ), center

  def get_config(self):
    config = {}
    base_config = super(ComputeInputStub, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return ( tf.shape(3,3), tf.shape(1,3) )

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
        # first, santity check that we know how to calc distances
        a = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2]
        ])
        b = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2]
        ])
        expected_dist = np.array([[1],
                                  [5],
                                  [6.232977 ]], dtype='float32')
        
        diff = a-b
        _, dist = tf.linalg.normalize( diff, axis=-1 )
        print( dist )
        np.testing.assert_array_almost_equal( dist, expected_dist )

        dist2 = tf.linalg.norm( diff, axis=-1, keepdims=True )
        np.testing.assert_array_almost_equal( dist2, expected_dist )
        
        ggp = a
        gp = b
        p = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2]
        ])

        layer = ComputeInputStub()
        output = layer([ ggp, gp, p ])

        # v
        np.testing.assert_array_almost_equal(
            output[1], 
            np.array([[20.  , 12.  , 12.  ],
                       [22.  ,  2.  , 42.  ],
                       [ 2.21,  2.2 , -2.2 ]], dtype='float32')
        )

        # M
        np.testing.assert_array_almost_equal(
            output[0], 
            np.array([[[ 0.78933704,  0.3803963 , -0.48191875],
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

class XRotMatRadians(Layer):
  def __init__(self,
               **kwargs):
    super(XRotMatRadians, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    angle_sin = tf.math.sin( inputs )
    angle_cos = tf.math.cos( inputs )
    
    ones = tf.ones_like( inputs )
    zeros = tf.zeros_like( inputs )

    angle_sin_neg = tf.math.negative( angle_sin );

    # doing columns instead of rows, even though the code reads like rows
    col_x = tf.concat( [ones,zeros,zeros], axis=-1 )
    col_y = tf.concat( [zeros,angle_cos,angle_sin], axis=-1 )
    col_z = tf.concat( [zeros,angle_sin_neg,angle_cos], axis=-1 )

    col_x = tf.expand_dims(col_x, -1)
    col_y = tf.expand_dims(col_y, -1)
    col_z = tf.expand_dims(col_z, -1)
    
    return tf.concat( [col_x,col_y,col_z], axis=-1 )

  def get_config(self):
    config = {}
    base_config = super(XRotMatRadians, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(1)

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
        a = np.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])

        layer = cls()
        output = layer( a ).numpy()
        np.testing.assert_array_almost_equal(
            output, 
            np.array([[[ 1.       , 0.       , 0.        ],
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


class ZRotMatRadians(Layer):
  def __init__(self,
               **kwargs):
    super(ZRotMatRadians, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    angle_sin = tf.math.sin( inputs )
    angle_cos = tf.math.cos( inputs )
    
    ones = tf.ones_like( inputs )
    zeros = tf.zeros_like( inputs )

    angle_sin_neg = tf.math.negative( angle_sin );

    # doing columns instead of rows, even though the code reads like rows
    col_x = tf.concat( [angle_cos,angle_sin,zeros], axis=-1 )
    col_y = tf.concat( [angle_sin_neg,angle_cos,zeros], axis=-1 )
    col_z = tf.concat( [zeros,zeros,ones], axis=-1 )

    col_x = tf.expand_dims(col_x, -1)
    col_y = tf.expand_dims(col_y, -1)
    col_z = tf.expand_dims(col_z, -1)
    
    return tf.concat( [col_x,col_y,col_z], axis=-1 )

  def get_config(self):
    config = { }
    base_config = super(ZRotMatRadians, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(1)

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
        a = np.array([
            [ 1., ],
            [ 2., ],
            [ 1.7, ],
            [ 0., ],
            [ -2.23, ],
        ])

        layer = cls()
        output = layer( a ).numpy()
        np.testing.assert_array_almost_equal(
            output, 
            np.array([[[ 0.5403023 , -0.84147096,  0.        ],
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


class TranslateFromStubB(Layer):
  def __init__(self,
               **kwargs):
    super(TranslateFromStubB, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    StubB = inputs[0]
    ic_in = inputs[1]

    col_x_indices = tf.constant([0,])
    col_x = tf.gather( StubB, col_x_indices, axis=2)
    col_x = tf.squeeze( col_x, axis=2 )

    d_indices = tf.constant([2,])
    d = tf.gather( ic_in, d_indices, axis=1 )
    
    vec = tf.math.multiply( d, col_x )
    return vec

  def get_config(self):
    config = {}
    base_config = super(TranslateFromStubB, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(1)

  @classmethod
  def debug(cls):
      #pylint: disable = E, W, R, C
      # Tested in GetAtomXYZ
      print( "passed {} tests".format( cls.__name__ ) )
      pass

class GetAtomXYZ(Layer):
  def __init__(self,
               **kwargs):
    super(GetAtomXYZ, self).__init__(**kwargs)
    self.comp_input_stub = ComputeInputStub()
    self.x_rot = XRotMatRadians()
    self.z_rot = ZRotMatRadians()
    self.trans = TranslateFromStubB()

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
    p_in = inputs[0]
    gp_in = inputs[1]
    ggp_in = inputs[2]
    ic_in = inputs[3]

    parent_M, parent_V = self.comp_input_stub([ggp_in,gp_in,p_in])

    phi = tf.gather( ic_in, tf.constant([0,]), axis=1 )
    x_phi = self.x_rot( phi )

    theta = tf.gather( ic_in, tf.constant([1,]), axis=1 )
    z_theta = self.z_rot( theta )

    StubA = tf.linalg.matmul( parent_M, x_phi )
    StubB = tf.linalg.matmul( StubA, z_theta )
    dv = self.trans([StubB,ic_in])
    v = parent_V + dv

    return v

  def get_config(self):
    config = {
    }
    base_config = super(GetAtomXYZ, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return tf.shape(3)

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
        ggp = np.array([
            [0., 1., 0.],
            [2., 4., 1.,],
            [2.6, -3.0, 0.2],
            [2.6, -3.0, 0.3],
        ])
        gp = np.array([
            [0., 1., 1.],
            [2., 0., 4.,],
            [2.1, 3.2, -0.2],
            [2.1, 3.2, -0.3]
        ])
        p = np.array([
            [20., 12., 12.],
            [22., 2., 42.,],
            [2.21, 2.2, -2.2],
            [2.21, 2.2, -2.3]
        ])
        ic = np.array([
            [0.3, -0.2, 0.1],
            [2.3, 0.2, 10],
            [3.4, 1.2, -0.2],
            [3.4, 1.2, 0.2]
        ])

        layer = cls()
        output = layer([ p, gp, ggp, ic ]).numpy()

        np.testing.assert_array_almost_equal(
            output, 
            np.array([[20.07297  , 12.033433 , 12.059646 ],
                       [24.958748 ,  1.6791693, 51.54688  ],
                       [ 2.1709127,  2.0679247, -2.0549886],
                       [ 2.2489984,  2.3320887, -2.445023 ]],
                      dtype='float32')
        )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class StubFromCoords(Layer):
  def __init__(self,
               **kwargs):
    super(StubFromCoords, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      N  = inputs[0]
      CA = inputs[1]
      C  = inputs[2]

      #matching vocab from Stub.cc
      a = N
      b = CA
      c = C
      center = CA

      e1 = a-b;
      e1, _ = tf.linalg.normalize( e1, axis=-1 )

      e3 = tf.linalg.cross( e1, c-b )
      e3, _ = tf.linalg.normalize( e3, axis=-1 )

      e2 = tf.linalg.cross( e3, e1 )

      e1 = tf.expand_dims(e1, -1)
      e2 = tf.expand_dims(e2, -1)
      e3 = tf.expand_dims(e3, -1)
      
      return tf.concat( [e1,e2,e3], axis=-1 ), center
      #return e2

  def get_config(self):
    config = { }
    base_config = super(StubFromCoords, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @tf_utils.shape_type_conversion
  def compute_output_shape(self, input_shape):
    return ( tf.shape(3,3), tf.shape(1,3) )

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
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

        o0 = np.array([[[ 0.        ,  0.87621593,  0.48191875],
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

        o1 = np.array([[ 0. ,  1. ,  1. ],
                        [ 2. ,  0. ,  4. ],
                        [ 2.1,  3.2, -0.2],
                        [ 2.1,  3.2, -0.3]], dtype='float32')

        layer = cls()
        output = [l.numpy() for l in layer([ N, CA, C ])]

        np.testing.assert_array_almost_equal( output[0], o0 )
        np.testing.assert_array_almost_equal( output[1], o1 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class HTValsFromCoords(Layer):
  def __init__(self, **kwargs):
    self.l = StubFromCoords()
    super(HTValsFromCoords, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      Ni  = inputs[0]
      CAi = inputs[1]
      Ci  = inputs[2]
      Nf  = inputs[3]
      CAf = inputs[4]
      Cf  = inputs[5]

      Mi, Vi = self.l([Ni,CAi,Ci])
      Mf, Vf = self.l([Nf,CAf,Cf])

      Mi_t = tf.transpose( Mi, perm=[0,2,1] )

      dV = Vf - Vi      
      dV = tf.expand_dims(dV, -1)
      
      T = tf.linalg.matmul( Mi_t, dV )
      R = tf.linalg.matmul( Mi_t, Mf )

      CAdist = tf.linalg.norm( CAf-CAi, axis=-1, keepdims=True )

      Tflat = Reshape(target_shape=[3])(T)
      Rflat = Reshape(target_shape=[9])(R)
      TRflat = tf.concat([Tflat,Rflat,CAdist], axis=-1)
      return TRflat

  def get_config(self):
    config = {}
    base_config = super(HTValsFromCoords, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
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
        o = np.array([[-1.1000000e+01,  2.2825424e+01,  0.0000000e+00,  4.3413538e-01,
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
        output = layer([ p0, p1, p2, p1, p2, p3 ]).numpy()
        #print( output )

        np.testing.assert_array_almost_equal( output, o )

        print( "passed {} tests".format( cls.__name__ ) )
        pass
    

class BuildHTs(Layer):
  def __init__(self,
               **kwargs):
    super(BuildHTs, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      A = inputs[0]

      N  = inputs[1]
      CA = inputs[2]
      C  = inputs[3]

      index_sources = A[:, 0]  # Nodes sending the message
      index_targets = A[:, 1]  # Nodes receiving the message

      Nsource = tf.gather( N, index_sources, axis=-2 )
      CAsource = tf.gather( CA, index_sources, axis=-2 )
      Csource = tf.gather( C, index_sources, axis=-2 )
      
      Ntarget = tf.gather( N, index_targets, axis=-2 )
      CAtarget = tf.gather( CA, index_targets, axis=-2 )
      Ctarget = tf.gather( C, index_targets, axis=-2 )

      return HTValsFromCoords()([Nsource,CAsource,Csource, Ntarget,CAtarget,Ctarget])

  def get_config(self):
    config = { }
    base_config = super(BuildHTs, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @classmethod
  def debug(cls):
        #pylint: disable = E, W, R, C
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
        output = layer([ A, N, CA, C ]).numpy()

        layer2 = HTValsFromCoords()
        o2 = layer2([ Ni, CAi, Ci, Nf, CAf, Cf ]).numpy()

        np.testing.assert_array_almost_equal( output, o2 )

        print( "passed {} tests".format( cls.__name__ ) )
        pass

class SameSide(Layer):
  def __init__(self,
               **kwargs):
    super(SameSide, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      A = inputs[0]
      NodeLabels  = inputs[1]

      index_sources = A.indices[:, 0]  # Nodes sending the message
      index_targets = A.indices[:, 1]  # Nodes receiving the message

      Lsource = tf.gather( NodeLabels, index_sources, axis=-2 )
      Ltarget = tf.gather( NodeLabels, index_targets, axis=-2 )

      L = tf.math.abs( Lsource - Ltarget )

      return L

  def get_config(self):
    config = { }
    base_config = super(SameSide, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

  @classmethod
  def debug(cls):
      #pylint: disable = E, W, R, C
      # TODO
      print( "passed {} tests".format( cls.__name__ ) )
      pass

class BuildIPAAtoms(Layer):
  def __init__(self,
               **kwargs):
    super(BuildIPAAtoms, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      N = inputs[0]
      CA = inputs[1]
      C = inputs[2]

      c1_icoor = tf.constant( [-1.57, 1.06, 1.01], shape=(1,3), dtype="float32" )
      c1_icoor = tf.repeat( c1_icoor, tf.shape(N)[0], axis=0 )
      child1 = GetAtomXYZ()([N, CA, C, c1_icoor ])

      c2_icoor = tf.constant( [0.52, 1.06, 1.01], shape=(1,3), dtype="float32" )
      c2_icoor = tf.repeat( c2_icoor, tf.shape(N)[0], axis=0 )
      child2 = GetAtomXYZ()([N, CA, C, c2_icoor ])

      c3_icoor = tf.constant( [2.61, 1.06, 1.01], shape=(1,3), dtype="float32" )
      c3_icoor = tf.repeat( c3_icoor, tf.shape(N)[0], axis=0 )
      child3 = GetAtomXYZ()([N, CA, C, c3_icoor ])

      # N is the parent
      return child1, child2, child3
            
  def get_config(self):
    config = { }
    base_config = super(BuildIPAAtoms, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

class IPADists(Layer):
  def __init__(self,
               **kwargs):
    super(IPADists, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      C = inputs[0]
      Nchild1 = inputs[1]
      Nchild2 = inputs[2]
      Nchild3 = inputs[3]
      A = inputs[4]
      N = inputs[5]

      index_sources = A.indices[:, 0]  # Nodes sending the message
      index_targets = A.indices[:, 1]  # Nodes receiving the message

      Cs = tf.gather( C, index_sources, axis=0 )
      target_children1 = tf.gather( Nchild1, index_targets, axis=0 )
      target_children2 = tf.gather( Nchild2, index_targets, axis=0 )
      target_children3 = tf.gather( Nchild3, index_targets, axis=0 )
      Ns = tf.gather( C, index_targets, axis=0 )
      
      dist1 = tf.linalg.norm( Cs-target_children1, axis=-1, keepdims=True )
      dist2 = tf.linalg.norm( Cs-target_children2, axis=-1, keepdims=True )
      dist3 = tf.linalg.norm( Cs-target_children3, axis=-1, keepdims=True )
      dist4 = tf.linalg.norm( Cs-Ns, axis=-1, keepdims=True )

      results = tf.concat([dist1,dist2,dist3,dist4], axis=-1)

      return results
      
  def get_config(self):
    config = { }
    base_config = super(IPADists, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))

class ExtractIPADists(Layer):
  def __init__(self,
               **kwargs):
    super(ExtractIPADists, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      N = inputs[0]
      CA = inputs[1]
      C = inputs[2]    
      A = inputs[3]
      child1, child2, child3 = BuildIPAAtoms()([N,CA,C])
      results = IPADists()([C, child1, child2, child3, A, N])
      return results
    
  def get_config(self):
    config = { }
    base_config = super(ExtractIPADists, self).get_config()
    return dict(list(base_config.items()) + list(config.items()))
  
class ExtractPeptideBondMeasurements(Layer):
  def __init__(self,
               **kwargs):
    super(ExtractPeptideBondMeasurements, self).__init__(**kwargs)

  @tf_utils.shape_type_conversion
  def build(self, input_shape):
    self.built = True

  def call(self, inputs):
      N = inputs[0]
      CA = inputs[1]
      C = inputs[2]    
      A = inputs[3]

      index_sources = A.indices[:, 0]  # Nodes sending the message
      index_targets = A.indices[:, 1]  # Nodes receiving the message
      Nsrc = tf.gather( N, index_sources, axis=0 )
      CAsrc = tf.gather( CA, index_sources, axis=0 )
      Csrc = tf.gather( C, index_sources, axis=0 )
      Ctgt = tf.gather( C, index_targets, axis=0 )
      CAtgt = tf.gather( CA, index_targets, axis=0 )
      Ntgt = tf.gather( N, index_targets, axis=0 )

      dist = tf.linalg.norm( Ctgt-Nsrc, axis=-1, keepdims=True )
      angle = Angle()([Ctgt,Nsrc,CAsrc])
      dihedral = Dihedral(apply_atan2=False)(
        [Ctgt,Nsrc,CAsrc,Csrc])

      angle2 = Angle()([CAtgt,Ctgt,Nsrc])
      dihedral2 = Dihedral(apply_atan2=False)(
        [CAtgt,Ctgt,Nsrc,CAsrc])

      dihedral3 = Dihedral(apply_atan2=False)(
        [Ntgt,CAtgt,Ctgt,Nsrc])
      
      return tf.concat( [dist,angle,dihedral,angle2,dihedral2,dihedral3], axis=-1 )
    
    
if __name__ == '__main__':
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
  SameSide.debug()
```