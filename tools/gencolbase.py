
# base colors:
# red, blue, teal, pink, brown, magenta, darkgray, violet, orange, olive
bcol = [[1,0,0],
      [0,0,1],
      [0,.5,.5],
      [1,.75,.75],
      [.75,.5,.25],
      [1,0,1],
      [.25,.25,.25],
      [.5,0,.5],
      [1,.5,0],
      [.5,.5,0]]

for i in range(0, 10):
    colr = int(255 * bcol[i][0])
    colg = int(255 * bcol[i][1])
    colb = int(255 * bcol[i][2])
    colbr = (colr + 255 * 2) / 3
    colbg = (colg + 255 * 2) / 3
    colbb = (colb + 255 * 2) / 3
    print('.mot{0} {{\n color: #{1:X}{2:X}{3:X};\n\
    background: #{4:X}{5:X}{6:X};\n\
    font-weight:bold;\n\
}}'.format(i+1, colr, colg, colb, colbr, colbg, colbb))
    if (i + 1) % 2:
        print('table tr:nth-child(even) td.mot{0} {{\n\
    background: #{1:X}{2:X}{3:X};\n\
}}'.format(i+1, colr, colg, colb))
    else:
        print('td.mot{0} {{\n\
    background: #{1:X}{2:X}{3:X};\n\
}}'.format(i+1, colr, colg, colb))

