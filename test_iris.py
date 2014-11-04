import iris
data_path = '/home/nicholat/project/pacemaker/um_output/'
u = iris.load(data_path + 'uagdsa@da040a1','eastward_wind')

# uagdsa@da040a1
# uagdsa@pa040ag
# uagdsa@pa040ar
# uagdsa@pa040dc
# uagdsa@pa040fb
# uagdsa@pa040ja
# uagdsa@pa040jl
# uagdsa@pa040jn
# uagdsa@pa040mr
# uagdsa@pa040my
# uagdsa@pa040nv
# uagdsa@pa040ot
# uagdsa@pa040sp
iris.io.save(u,'utest.nc')

