import iris
import time

print time.ctime()
print 'load all_var'
all_var = iris.load('um_output/uagdsa@pa*')
print 'finished loading all'
print time.ctime()
print all_var[44]
print time.ctime()
print 'load some_var'
some_var = iris.load('um_output/uagdsa@pa*')[44]
print 'finished loading some'
print time.ctime()
print some_var
print time.ctime()
print 'load some_var'
some_var = iris.load('um_output/uagdsa@pa*','x_wind')[0]
print 'finished loading some'
print time.ctime()
print some_var
