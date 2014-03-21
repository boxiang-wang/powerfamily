      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1223 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1224
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1225
      integer jd(*),ia(nx),nin(nlam)                                       1226
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11831                                    1230
      jerr=10000                                                           1230
      return                                                               1230
11831 continue                                                             1231
      allocate(ww(1:no),stat=jerr)                                         1232
      allocate(ju(1:ni),stat=ierr)                                         1232
      jerr=jerr+ierr                                                       1233
      allocate(vq(1:ni),stat=ierr)                                         1233
      jerr=jerr+ierr                                                       1234
      allocate(xm(1:ni),stat=ierr)                                         1234
      jerr=jerr+ierr                                                       1235
      if(isd .le. 0)goto 11851                                             1235
      allocate(xs(1:ni),stat=ierr)                                         1235
      jerr=jerr+ierr                                                       1235
11851 continue                                                             1236
      if(jerr.ne.0) return                                                 1237
      call chkvars(no,ni,x,ju)                                             1238
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1239
      if(maxval(ju) .gt. 0)goto 11871                                      1239
      jerr=7777                                                            1239
      return                                                               1239
11871 continue                                                             1240
      vq=max(0.0,vp)                                                       1240
      vq=vq*ni/sum(vq)                                                     1241
11880 do 11881 i=1,no                                                      1241
      ww(i)=sum(y(i,:))                                                    1241
      y(i,:)=y(i,:)/ww(i)                                                  1241
11881 continue                                                             1241
11882 continue                                                             1241
      sw=sum(ww)                                                           1241
      ww=ww/sw                                                             1242
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1243
      if(nc .ne. 1)goto 11901                                              1244
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1246 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11911                                                           1247
11901 continue                                                             1248
11911 continue                                                             1251
11891 continue                                                             1251
      if(jerr.gt.0) return                                                 1251
      dev0=2.0*sw*dev0                                                     1252
11920 do 11921 k=1,lmu                                                     1252
      nk=nin(k)                                                            1253
11930 do 11931 ic=1,nc                                                     1253
      if(isd .le. 0)goto 11951                                             1253
11960 do 11961 l=1,nk                                                      1253
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1253
11961 continue                                                             1253
11962 continue                                                             1253
11951 continue                                                             1254
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1255
11931 continue                                                             1256
11932 continue                                                             1256
11921 continue                                                             1257
11922 continue                                                             1257
      deallocate(ww,ju,vq,xm)                                              1257
      if(isd.gt.0) deallocate(xs)                                          1258
      return                                                               1259
      end                                                                  1260
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1261
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1261
      integer ju(ni)                                                       1262
11970 do 11971 j=1,ni                                                      1262
      if(ju(j).eq.0)goto 11971                                             1263
      xm(j)=dot_product(w,x(:,j))                                          1263
      x(:,j)=x(:,j)-xm(j)                                                  1264
      if(isd .le. 0)goto 11991                                             1264
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1264
      x(:,j)=x(:,j)/xs(j)                                                  1264
11991 continue                                                             1265
11971 continue                                                             1266
11972 continue                                                             1266
      return                                                               1267
      end                                                                  1268
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1270 
     *,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1272 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1273
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1274
      integer ju(ni),m(nx),kin(nlam)                                       1275
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1280
      allocate(xv(1:ni),stat=ierr)                                         1280
      jerr=jerr+ierr                                                       1281
      allocate(ga(1:ni),stat=ierr)                                         1281
      jerr=jerr+ierr                                                       1282
      allocate(bs(0:ni),stat=ierr)                                         1282
      jerr=jerr+ierr                                                       1283
      allocate(mm(1:ni),stat=ierr)                                         1283
      jerr=jerr+ierr                                                       1284
      allocate(ixx(1:ni),stat=ierr)                                        1284
      jerr=jerr+ierr                                                       1285
      allocate(r(1:no),stat=ierr)                                          1285
      jerr=jerr+ierr                                                       1286
      allocate(v(1:no),stat=ierr)                                          1286
      jerr=jerr+ierr                                                       1287
      allocate(q(1:no),stat=ierr)                                          1287
      jerr=jerr+ierr                                                       1288
      if(jerr.ne.0) return                                                 1289
      fmax=log(1.0/pmin-1.0)                                               1289
      fmin=-fmax                                                           1289
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1290
      bta=parm                                                             1290
      q0=dot_product(w,y)                                                  1291
      if(q0 .gt. pmin)goto 12011                                           1291
      jerr=8001                                                            1291
      return                                                               1291
12011 continue                                                             1292
      if(q0 .lt. 1.0-pmin)goto 12031                                       1292
      jerr=9001                                                            1292
      return                                                               1292
12031 continue                                                             1293
      ixx=1                                                                1293
      al=0.0                                                               1293
      bz=log(q0/(1.0-q0))                                                  1294
      if(nonzero(no,g) .ne. 0)goto 12051                                   1294
      vi=q0*(1.0-q0)                                                       1294
      b(0)=bz                                                              1294
      v=vi*w                                                               1295
      r=w*(y-q0)                                                           1295
      q=q0                                                                 1295
      xmz=vi                                                               1295
      dev1=-(bz*q0+log(1.0-q0))                                            1296
      goto 12061                                                           1297
12051 continue                                                             1297
      b(0)=azero(no,y,g,w,jerr)                                            1297
      if(jerr.ne.0) return                                                 1298
      q=1.0/(1.0+exp(-b(0)-g))                                             1298
      v=w*q*(1.0-q)                                                        1298
      r=w*(y-q)                                                            1298
      xmz=sum(v)                                                           1299
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1300
12061 continue                                                             1301
12041 continue                                                             1301
      if(kopt .le. 0)goto 12081                                            1302
      if(isd .le. 0)goto 12101                                             1302
      xv=0.25                                                              1302
      goto 12111                                                           1303
12101 continue                                                             1303
12120 do 12121 j=1,ni                                                      1303
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1303
12121 continue                                                             1303
12122 continue                                                             1303
12111 continue                                                             1304
12091 continue                                                             1304
12081 continue                                                             1305
      dev0=dev1                                                            1306
12130 do 12131 i=1,no                                                      1306
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1307
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1308
12131 continue                                                             1309
12132 continue                                                             1309
      if(flmin .ge. 1.0)goto 12151                                         1309
      eqs=max(eps,flmin)                                                   1309
      alf=eqs**(1.0/(nlam-1))                                              1309
12151 continue                                                             1310
      m=0                                                                  1310
      mm=0                                                                 1310
      nlp=0                                                                1310
      nin=nlp                                                              1310
      mnl=min(mnlam,nlam)                                                  1310
      bs=0.0                                                               1310
      b(1:ni)=0.0                                                          1311
      shr=shri*dev0                                                        1312
12160 do 12161 j=1,ni                                                      1312
      if(ju(j).eq.0)goto 12161                                             1312
      ga(j)=abs(dot_product(r,x(:,j)))                                     1312
12161 continue                                                             1313
12162 continue                                                             1313
12170 do 12171 ilm=1,nlam                                                  1313
      al0=al                                                               1314
      if(flmin .lt. 1.0)goto 12191                                         1314
      al=ulam(ilm)                                                         1314
      goto 12181                                                           1315
12191 if(ilm .le. 2)goto 12201                                             1315
      al=al*alf                                                            1315
      goto 12181                                                           1316
12201 if(ilm .ne. 1)goto 12211                                             1316
      al=big                                                               1316
      goto 12221                                                           1317
12211 continue                                                             1317
      al0=0.0                                                              1318
12230 do 12231 j=1,ni                                                      1318
      if(ju(j).eq.0)goto 12231                                             1318
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1318
12231 continue                                                             1319
12232 continue                                                             1319
      al=alf*al0                                                           1320
12221 continue                                                             1321
12181 continue                                                             1321
      tlam=(2.0*al-al0)                                                       
! 12240 do 12241 k=1,ni                                                      1322
!       if(ixx(k).eq.1)goto 12241                                            1322
!       if(ju(k).eq.0)goto 12241                                             1323
!       if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1324
! 12241 continue                                                             1325
12242 continue                                                             1325
10680 continue                                                             1326
12250 continue                                                             1326
12251 continue                                                             1326
      bs(0)=b(0)                                                           1326
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1327
      if(kopt .ne. 0)goto 12271                                            1328
12280 do 12281 j=1,ni                                                      1328
!       if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1328
      xv(j)=dot_product(v,x(:,j)**2)                                       1328
12281 continue                                                             1329
12282 continue                                                             1329
12271 continue                                                             1330
12290 continue                                                             1330
12291 continue                                                             1330
      nlp=nlp+1                                                            1330
      dlx=0.0                                                              1331
12300 do 12301 k=1,ni                                                      1331
!       if(ixx(k).eq.0)goto 12301                                            1332
      bk=b(k)                                                              1332
      gk=dot_product(r,x(:,k))                                             1333
      u=gk+xv(k)*b(k)                                                      1333
      au=abs(u)-vp(k)*al
      if(au .gt. 0.0)goto 12321                                            1334
      b(k)=0.0                                                             1334
      goto 12331                                                           1335
12321 continue                                                             1335
      b(k)=sign(au,u)/(xv(k)+vp(k)*bta)                                    1335
12331 continue                                                             1336
12311 continue                                                             1336
      d=b(k)-bk                                                            1336
      if(abs(d).le.0.0)goto 12301                                          1336
      dlx=max(dlx,xv(k)*d**2)                                              1337
      r=r-d*v*x(:,k)                                                       1338
      if(mm(k) .ne. 0)goto 12351                                           1338
      nin=nin+1                                                            1338
      if(nin.gt.nx)goto 12302                                              1339
      mm(k)=nin                                                            1339
      m(nin)=k                                                             1340
12351 continue                                                             1341
12301 continue                                                             1342
12302 continue                                                             1342
      if(nin.gt.nx)goto 12292                                              1343
      d=sum(r)/xmz                                                         1344
      if(d .eq. 0.0)goto 12371                                             1344
      b(0)=b(0)+d                                                          1344
      dlx=max(dlx,xmz*d**2)                                                1344
      r=r-d*v                                                              1344
12371 continue                                                             1345
      if(dlx.lt.shr)goto 12292                                             1345
      if(nlp .le. maxit)goto 12391                                         1345
      jerr=-ilm                                                            1345
      return                                                               1345
12391 continue                                                             1346
12400 continue                                                             1346
12401 continue                                                             1346
      nlp=nlp+1                                                            1346
      dlx=0.0                                                              1347
12410 do 12411 l=1,nin                                                     1347
      k=m(l)                                                               1347
      bk=b(k)                                                              1348
      gk=dot_product(r,x(:,k))                                             1349
      u=gk+xv(k)*b(k)                                                      1349
      au=abs(u)-vp(k)*al                                                 1350
      if(au .gt. 0.0)goto 12431                                            1350
      b(k)=0.0                                                             1350
      goto 12441                                                           1351
12431 continue                                                             1351
      b(k)=sign(au,u)/(xv(k)+vp(k)*bta)                                    1351
12441 continue                                                             1352
12421 continue                                                             1352
      d=b(k)-bk                                                            1352
      if(abs(d).le.0.0)goto 12411                                          1352
      dlx=max(dlx,xv(k)*d**2)                                              1353
      r=r-d*v*x(:,k)                                                       1354
12411 continue                                                             1355
12412 continue                                                             1355
      d=sum(r)/xmz                                                         1356
      if(d .eq. 0.0)goto 12461                                             1356
      b(0)=b(0)+d                                                          1356
      dlx=max(dlx,xmz*d**2)                                                1356
      r=r-d*v                                                              1356
12461 continue                                                             1357
      if(dlx.lt.shr)goto 12402                                             1357
      if(nlp .le. maxit)goto 12481                                         1357
      jerr=-ilm                                                            1357
      return                                                               1357
12481 continue                                                             1358
      goto 12401                                                           1359
12402 continue                                                             1359
      goto 12291                                                           1360
12292 continue                                                             1360
      if(nin.gt.nx)goto 12252                                              1361
12490 do 12491 i=1,no                                                      1361
      fi=b(0)+g(i)                                                         1362
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1363
      if(fi .ge. fmin)goto 12511                                           1363
      q(i)=0.0                                                             1363
      goto 12501                                                           1363
12511 if(fi .le. fmax)goto 12521                                           1363
      q(i)=1.0                                                             1363
      goto 12531                                                           1364
12521 continue                                                             1364
      q(i)=1.0/(1.0+exp(-fi))                                              1364
12531 continue                                                             1365
12501 continue                                                             1365
12491 continue                                                             1366
12492 continue                                                             1366
      v=w*q*(1.0-q)                                                        1366
      xmz=sum(v)                                                           1366
      if(xmz.le.vmin)goto 12252                                            1366
      r=w*(y-q)                                                            1367
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12551                           1367
      ix=0                                                                 1368
12560 do 12561 j=1,nin                                                     1368
      k=m(j)                                                               1369
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12561                           1369
      ix=1                                                                 1369
      goto 12562                                                           1370
12561 continue                                                             1371
12562 continue                                                             1371
      if(ix .ne. 0)goto 12581                                              1372
12590 do 12591 k=1,ni                                                      1372
!       if(ixx(k).eq.1)goto 12591                                            1372
!       if(ju(k).eq.0)goto 12591                                             1373
      ga(k)=abs(dot_product(r,x(:,k)))                                     1374
!       if(ga(k) .le. al*vp(k))goto 12611                                   1374
!       ixx(k)=1                                                             1374
!       ix=1                                                                 1374
12611 continue                                                             1375
12591 continue                                                             1376
12592 continue                                                             1376
!       if(ix.eq.1) go to 10680                                              1377
      goto 12252                                                           1378
12581 continue                                                             1379
12551 continue                                                             1380
      goto 12251                                                           1381
12252 continue                                                             1381
      if(nin .le. nx)goto 12631                                            1381
      jerr=-10000-ilm                                                      1381
      goto 12172                                                           1381
12631 continue                                                             1382
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1382
      kin(ilm)=nin                                                         1383
      a0(ilm)=b(0)                                                         1383
      alm(ilm)=al                                                          1383
      lmu=ilm                                                              1384
      devi=dev2(no,w,y,q,pmin)                                             1385
      dev(ilm)=(dev1-devi)/dev0                                            1385
      if(xmz.le.vmin)goto 12172                                            1386
      if(ilm.lt.mnl)goto 12171                                             1386
      if(flmin.ge.1.0)goto 12171                                           1387
      me=0                                                                 1387
12640 do 12641 j=1,nin                                                     1387
      if(a(j,ilm).ne.0.0) me=me+1                                          1387
12641 continue                                                             1387
12642 continue                                                             1387
      if(me.gt.ne)goto 12172                                               1388
      if(dev(ilm).gt.devmax)goto 12172                                     1388
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12172                             1389
12171 continue                                                             1390
12172 continue                                                             1390
      g=log(q/(1.0-q))                                                     1391
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1392
      return                                                               1393
      end                                                                  1394
      function dev2(n,w,y,p,pmin)                                          1395
      real w(n),y(n),p(n)                                                  1396
      pmax=1.0-pmin                                                        1396
      s=0.0                                                                1397
12650 do 12651 i=1,n                                                       1397
      pi=min(max(pmin,p(i)),pmax)                                          1398
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1399
12651 continue                                                             1400
12652 continue                                                             1400
      dev2=s                                                               1401
      return                                                               1402
      end                                                                  1403
      function azero(n,y,g,q,jerr)                                         1404
      parameter(eps=1.0e-7)                                                1405
      real y(n),g(n),q(n)                                                  1406
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1410
      allocate(p(1:n),stat=ierr)                                           1410
      jerr=jerr+ierr                                                       1411
      allocate(w(1:n),stat=ierr)                                           1411
      jerr=jerr+ierr                                                       1412
      if(jerr.ne.0) return                                                 1413
      az=0.0                                                               1413
      e=exp(-g)                                                            1413
      qy=dot_product(q,y)                                                  1413
      p=1.0/(1.0+e)                                                        1414
12660 continue                                                             1414
12661 continue                                                             1414
      w=q*p*(1.0-p)                                                        1415
      d=(qy-dot_product(q,p))/sum(w)                                       1415
      az=az+d                                                              1415
      if(abs(d).lt.eps)goto 12662                                          1416
      ea0=exp(-az)                                                         1416
      p=1.0/(1.0+ea0*e)                                                    1417
      goto 12661                                                           1418
12662 continue                                                             1418
      azero=az                                                             1419
      deallocate(e,p,w)                                                    1420
      return                                                               1421
      end                                                                  1422
      subroutine chkvars(no,ni,x,ju)                                        888
      real x(no,ni)                                                         888
      integer ju(ni)                                                        889
10860 do 10861 j=1,ni                                                       889
      ju(j)=0                                                               889
      t=x(1,j)                                                              890
10870 do 10871 i=2,no                                                       890
      if(x(i,j).eq.t)goto 10871                                             890
      ju(j)=1                                                               890
      goto 10872                                                            890
10871 continue                                                              891
10872 continue                                                              891
10861 continue                                                              892
10862 continue                                                              892
      return                                                                893
      end                                                                   894
      function nonzero(n,v)                                                2494
      real v(n)                                                            2495
      nonzero=0                                                            2495
17010 do 17011 i=1,n                                                       2495
      if(v(i) .eq. 0.0)goto 17031                                          2495
      nonzero=1                                                            2495
      return                                                               2495
17031 continue                                                             2495
17011 continue                                                             2496
17012 continue                                                             2496
      return                                                               2497
      end                                                                  2498