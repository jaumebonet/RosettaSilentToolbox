# @Author: Jaume Bonet <bonet>
# @Date:   30-Jun-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: utils.py
# @Last modified by:   bonet
# @Last modified time: 30-Jun-2017


def selection2list( strrange ):
    r = strrange.split(",")
    o = []
    for x in r:
        if "-" not in x:
            o.append(int(x))
        else:
            xx = x.split("-")
            for i in range(int(xx[0]),  int(xx[1]) + 1):
                o.append(i)
    return o
