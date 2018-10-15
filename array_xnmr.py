data_dir = '/home/manning/nmr/nmr_data/'
fname = 'kimberley/20150529-4/params'

f = open(data_dir+fname,'r')
x = f.readlines()
f.close()
found=-1
for i in range(len(x)):
    if x[i][0] == ';':
        found = i
        break

print(found)
fout = open('params_out','w')
fout.writelines(x[0:found])

pw1 = np.array([0.050, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000])*1e-3

for i in range(1,4*len(pw1)):
    fout.write(';\n')
    if i%4 == 0:
        fout.write(('pw1 = %1.3fe-3\n' % (pw1[int(i/4)]*1000)))
    fout.write(('exp_index = %i\n' % i))
fout.close()


