import os, sys

def remesh( obj, size, feat, iters=10 ):
    # create the output path
    if not os.path.exists('output'):
        os.mkdir('output')
    # generate the remesher command line
    cmd = 'bin/remesher -input meshes/%s.obj -output output/%s_remesh.obj -size %f -feat %f -iters %d' % (obj,obj,size,feat,iters)
    print cmd
    os.system(cmd)

#remesh( 'u_joint_link',            0.02,  80.0, 20 )
#remesh( 'axis_end',                0.02,  30.0, 20 )
#remesh( 'bunny',                   0.01, 180.0, 20 )
#remesh('phone_mount',              0.02,  30.0, 40 )
#remesh( 'lever_arm',               0.01,  30.0, 20 )
#remesh( 'delta_arm_end_effector',  0.02,  30.0, 20 )
#remesh( 'chain_fixed',             0.03,  30.0, 20 )
#remesh('anti_backlash_nut',        0.02,  30.0, 20 )
#remesh( 'cap_thing',               0.02,  30.0, 20 )
remesh( 'delta_arm2',              0.015,  30.0, 20 )

sys.exit()

#remesh( 'axis_end',                0.01,  30.0, 5 )
#remesh( 'bearing_plate_608',       0.05,  30.0 )
#remesh( 'carriage',                0.05,  30.0 )
#remesh( 'delta_arm_base',          0.05,  30.0 )
#remesh( 'gear2',                   0.02,  30.0 )
#remesh('phone_mount',              0.01,  30.0 )
sys.exit()
remesh( 'shaft_support_8mm',       0.05,  30.0 )
remesh( 'u_joint_link',            0.05,  30.0 )
remesh( 'delta_robot_linear_axis', 0.05,  30.0 )

remesh('9x1_ratio_plate',          0.03,  30.0 )

