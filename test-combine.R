x <- BS.cancer.ex.fit[1:1000, 1]
y <- BS.cancer.ex.fit[501:1500, 2]
z <- BS.cancer.ex.fit[701:1700, 3:4]
xx <- BS.cancer.ex.fit.hdf5[1:1000, 1]
yy <- BS.cancer.ex.fit.hdf5[501:1500, 2]
zz <- BS.cancer.ex.fit.hdf5[701:1700, 3:4]

a <- combine(x, y)
aa <- combine(xx, yy)
b <- combine(x, y, z)
bb <- combine(xx, yy, zz)

a1 <- combineList(x, y)
aa1 <- combineList(xx, yy)
b1 <- combineList(x, y, z)
bb1 <- combineList(xx, yy, zz)

# TODO: a and a1 not identical because combine() preserves @trans even when
#       objects have different GRanges
all.equal(a, a1)

u <- BS.cancer.ex.fit[, 1]
v <- BS.cancer.ex.fit[, 2]
w <- BS.cancer.ex.fit[, 3:4]
uu <- BS.cancer.ex.fit.hdf5[, 1]
vv <- BS.cancer.ex.fit.hdf5[, 2]
ww <- BS.cancer.ex.fit.hdf5[, 3:4]

c <- combine(u, v)
cc <- combine(uu, vv)
d <- combine(u, v, w)
dd <- combine(uu, vv, ww)

c1 <- combineList(u, v)
cc1 <- combineList(uu, vv)
d1 <- combineList(u, v, w)
dd1 <- combineList(uu, vv, ww)

# TODO: FALSE; investigate
all.equal(c, c1)
