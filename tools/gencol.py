import random



for i in range(1, 11):
    colr = random.randint(128, 255)
    colg = random.randint(128, 255)
    colb = random.randint(128, 255)
    colbr = (colr + 255 * 2) / 3
    colbg = (colg + 255 * 2) / 3
    colbb = (colb + 255 * 2) / 3
    print('.mot{0} {{\n color: #{1:X}{2:X}{3:X};\n\
    background: #{4:X}{5:X}{6:X};\n\
    font-weight:bold;\n\
}}'.format(i, colr, colg, colb, colbr, colbg, colbb))
    if i % 2:
        print('table tr:nth-child(even) td.mot{0} {{\n\
    background: #{1:X}{2:X}{3:X};\n\
}}'.format(i, colbr, colbg, colbb))
    else:
        print('td.mot{0} {{\n\
    background: #{1:X}{2:X}{3:X};\n\
}}'.format(i, colbr, colbg, colbb))

