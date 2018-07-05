package net.arwix.astronomy2.ephemeris.moshier

/** heliocentric earth-moon barycenter
 *   polar coordinates re fixed J2000 mean equinox and ecliptic
 *   S. L. Moshier
 *   December, 1996
 *
 *              Residuals against JPL ephemeris DE404 (arc seconds)
 * First date in file = 1221000.50
 * Number of samples = 524290
 * Sampling interval = 3.0 days
 * Julian Years             Longitude          Latitude           Distance
 *                                                               1 = 725 km
 *                       Peak  RMS   Ave    Peak  RMS   Ave    Peak  RMS   Ave
 * -1369.0 to -1000.0:   0.07  0.02  0.00   0.02  0.00 -0.00   0.03  0.01 -0.00
 * -1000.0 to  -500.0:   0.05  0.01 -0.00   0.01  0.00 -0.00   0.03  0.01 -0.00
 *  -500.0 to     0.0:   0.04  0.01  0.00   0.01  0.00  0.00   0.03  0.01 -0.00
 *     0.0 to   500.0:   0.05  0.01 -0.00   0.01  0.00  0.00   0.03  0.01 -0.00
 *   500.0 to  1000.0:   0.06  0.01  0.00   0.01  0.00  0.00   0.03  0.01 -0.00
 *  1000.0 to  1500.0:   0.05  0.01 -0.00   0.01  0.00  0.00   0.02  0.01 -0.00
 *  1500.0 to  2000.0:   0.05  0.01  0.00   0.01  0.00  0.00   0.02  0.00  0.00
 *  2000.0 to  2500.0:   0.05  0.01 -0.00   0.01  0.00 -0.00   0.03  0.01 -0.00
 *  2500.0 to  2937.2:   0.07  0.02  0.00   0.02  0.00 -0.00   0.03  0.01 -0.00
 */
internal object EarthMoonBarycenterMoshierData : MoshierData() {

    override val tabl by lazy { doubleArrayOf(-242809.0, -178223.0, -6154.0, -6547.0, 15526.0, -79460.0, 66185.0, -19531.0, -12754.0, 4389.0, 3153.0, -1151.0, 768.0, 1750.0, -248.0, 657.0, -80.0, 0.0, -4.0, -29.0, -3020.0, 301.0, -360.0, 412.0, -1463.0, 2266.0, -41.0, 30.0, -39868.0, -14275.0, -25052.0, 1583.0, 15695.0, 10018.0, -113.0, -122.0, -243.0, 18.0, -33.0, 31.0, -134.0, -171.0, 243.0, -115.0, 18.0, 148.0, -120.0, -129.0, 19.0, -220.0, -30.0, 19.0, 8.0, 23.0, -162.0, -124.0, 189.0, -315.0, 73.0, 77.0, 32006.0, -11295.0, 11595.0, 5629.0, -838.0, 1728.0, 0.0, 4.0, 38.0, 15.0, 142.0, -228.0, 92.0, 32.0, -2274.0, -1500.0, -2277.0, 3143.0, 3204.0, 127.0, -20.0, -11.0, 5186.0, 1054.0, 996.0, 1769.0, -231.0, 163.0, -88.0, -19.0, -2.0, -145.0, -27.0, 48.0, -8.0, 421.0, -7.0, 148.0, -16.0, -2.0, -3964.0, 4259.0, -11192.0, -8385.0, 11513.0, -13415.0, 103.0, -43.0, -289.0, -79.0, -29.0, 163.0, -117.0, 559.0, -190.0, -15.0, 7108.0, 5345.0, 12933.0, -7709.0, 3485.0, -26023.0, 11.0, -5.0, 311.0, 78.0, 22.0, 76.0, 2846.0, -3922.0, 2329.0, 43.0, 34.0, 442.0, 3.0, -245.0, -5.0, -3.0, -17.0, 5.0, 318.0, 15963.0, 2520.0, 7115.0, 2548.0, -9836.0, -7063.0, 1950.0, -4471.0, -8326.0, 4964.0, -3101.0, 563.0, -80.0, -1444.0, -472.0, 8.0, -22.0, 1558.0, -88.0, 235.0, 359.0, 293.0, -16.0, 144.0, 209.0, -13.0, -7.0, 812.0, -744.0, 150.0, -740.0, -2569.0, -956.0, 69.0, -2475.0, 1009.0, -55.0, -1707.0, -2525.0, 1071.0, -1761.0, 550.0, 279.0, -14.0, 36.0, -10442.0, 3344.0, -6759.0, -21551.0, 24737.0, -434.0, 504.0, -385.0, 191.0, 96.0, -2760.0, -1068.0, 85.0, -2617.0, 1062.0, -43.0, 192.0, -16.0, 30.0, 42.0, -2446.0, 588.0, -1522.0, -2933.0, 1746.0, -1070.0, 511.0, -1401.0, 139.0, 729.0, -12.0, 29.0, -2618.0, -2076.0, 2079.0, -3711.0, -21.0, -2727.0, -80.0, -19.0, 113.0, 2420.0, 325.0, 1058.0, 379.0, -1478.0, 296.0, -251.0, -265.0, -409.0, -10.0, 20.0, 15.0, -15.0, 11.0, 143.0, -83.0, 19.0, 266.0, -17.0, 40.0, 59.0, 19.0, -105.0, 5.0, 48331.0, 21.0, -16.0, -97.0, -318.0, 158.0, -171.0, 456.0, -351.0, 168.0, 85.0, 12.0, -2.0, 20.0, -15.0, 15.0, 2.0, 385.0, -1125.0, 521.0, -23.0, -815.0, -2088.0, 1644.0, -1329.0, 7.0, 14.0, -582.0, 234.0, -67.0, -467.0, -167.0, -51.0, -684.0, -2302.0, 1315.0, -797.0, 6.0, -70.0, -118.0, -406.0, 67.0, -63.0, -4848.0, 3713.0, -8483.0, -8776.0, 13049.0, -9404.0, -23.0, 34.0, -12.0, 1.0, -24.0, -10.0, -21.0, 0.0, -1.0, 24.0, -3.0, 28.0, -3032.0, -2494.0, 2498.0, -4342.0, -6538.0, 1899.0, -4414.0, -13249.0, 15540.0, -292.0, -228.0, 176.0, -40.0, -161.0, -20.0, -36.0, -800.0, -172.0, -36.0, -208.0, -249.0, -374.0, -1410.0, -72118.0, -745.0, 213.0, -23.0, 196.0, -14.0, -2.0, -239.0, -341.0, 1015.0, -291.0, 33.0, -94.0, 90.0, -20431.0, 4.0, -39.0, 75.0, 216.0, -23.0, 41.0, 116.0, 24.0, 5.0, 26.0, -45.0, -4178.0, -9.0, -23.0, 12.0, 18.0, 68.0, -2.0, 36.0, -19.0, 42.0, -8.0, 6.0, -106.0, 4.0, -38.0, -73.0, 259.0, 107.0, -293.0, -12.0, -44.0, 37.0, 13.0, 73.0, -46.0, 17.0, 8.0, 5832.0, 1989.0, -1404.0, 4469.0, -1619.0, -743.0, -1796.0, -2206.0, 461.0, -291.0, 153.0, 1104.0, 19195.0, 652503.0, 5587.0, -5252787.0, 47.0, -17340051.0, -32.0, 68926621.0, 1054.0, -230.0, -1601.0, 356.0, -562.0, -998.0, 124.0, -446.0, -171.0, 66.0, 26.0, 60.0, -7.0, 0.0, -88.0, -43.0, 65.0, -400.0, 4.0, 183.0, -1014.0, 2043.0, -1076.0, 126.0, -41.0, -205.0, -127.0, -85.0, -15.0, 68.0, 0.0, 0.0, -320.0, 75.0, -42.0, 285.0, -303.0, 771.0, 616.0, 400.0, -470.0, 48.0, -76.0, -103.0, -190.0, -11.0, -139.0, -5.0, -48.0, -87.0, -22.0, -362.0, -271.0, 1233.0, -392.0, 353.0, -154.0, -71.0, -109.0, 112.0, 17.0, 8.0, 1.0, -17.0, -170.0, 623.0, -279.0, 21.0, 139.0, -151.0, -107.0, -55199.0, 588.0, -188.0, 397.0, 674.0, -406.0, 269.0, 166.0, -207.0, 585.0, 333.0, -386.0, 754.0, 29.0, -65.0, 35.0, 10.0, 63.0, 1291.0, 62.0, 8.0, 239.0, 1323.0, -1434.0, 53.0, 19.0, -1.0, 34.0, 82.0, -15.0, -16.0, 265.0, -338.0, -729.0, -207.0, 3.0, 17.0, 697.0, 399.0, 274.0, 760.0, -12.0, 2.0, -48.0, -9.0, 3.0, 64.0, 147.0, 36.0, 9.0, 46.0, 77.0, 144.0, -76.0, 65.0, 2329.0, 1763.0, 987.0, 5506.0, 66.0, -123.0, -41.0, -24.0, -12.0, 1.0, -19.0, 94.0, 19.0, 8.0, -1.0, -18.0, 142.0, 77.0, -78.0, 187.0, 6.0, 18.0, 607.0, 163.0, 17.0, 158.0, 27.0, -208.0, 154.0, 27317.0, 587.0, -143.0, 22.0, -153.0, 5.0, -34.0, 75.0, 330.0, 98.0, -128.0, -67.0, -6542.0, -115.0, -236.0, 217.0, -12.0, 10.0, -6.0, -250.0, 653.0, 1611.0, -209.0, 4.0, 1072.0, -129.0, 216.0, 402.0, 81.0, 117.0, 11.0, 0.0, 20.0, 24.0, -28.0, 245.0, 437.0, -16.0, 59.0, 527952.0, -74895.0, 169682.0, 177186.0, -376.0, -362869.0, -60.0, 719658.0, -151.0, -382.0, -22.0, -43.0, 5.0, -5.0, 14.0, 5.0, -9.0, 13.0, 83.0, 296.0, -369.0, -1.0, -14.0, -6.0, 42.0, 8.0, -31.0, 7.0, -354.0, 634.0, 1132.0, 243.0, -38.0, 42.0, -14.0, 68.0, -6.0, 31.0, -36.0, -13.0, 7.0, -2104.0, 16.0, 67.0, 9.0, -4.0, 174.0, 144.0, 58.0, 438.0, -15.0, 5.0, -16.0, 19.0, -135.0, 1642.0, -140.0, -11.0, -4.0, 27.0, 253.0, -382.0, -651.0, -221.0, 11.0, 1.0, -12.0, -6.0, 136.0, 23.0, -1.0, 43.0, 3.0, 38.0, -26.0, -5.0, 17864.0, -4244.0, 5704.0, 7754.0, -36.0, -7891.0, -3.0, 10418.0, 2.0, -844.0, -1.0, 126.0, -7.0, 32.0, -67.0, -5.0, 39.0, 10.0, 5.0, 52.0, -13.0, 159.0, -49.0, -21.0, 1.0, -394.0, 7.0, -15.0, -4.0, -245.0, 1.0, 172.0, -36.0, -3.0, 13.0, 5.0, 0.0, 1.0, -1.0, 0.0, 0.0, -202.0, -2.0, 19.0, -20.0, -2.0, 5.0, 3.0, 0.0, -110.0, -12.0, -1.0, 0.0, -62.0, 0.0, -36.0, 0.0, -22.0, -13.0, 3.0) }
    override val tabb by lazy { doubleArrayOf(-428091.0, -488399.0, 746850.0, 6.0, 210.0, -93.0, 32.0, 1.0, -365.0, 332.0, -105.0, 76.0, -7.0, 2.0, -8.0, 14.0, -1.0, 2.0, 0.0, 0.0, -65.0, 12.0, -17.0, 7.0, -1.0, 1.0, 0.0, 0.0, -15.0, 65.0, -4.0, 26.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, -1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -1.0, 0.0, -30.0, 28.0, -6.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, -16.0, 20.0, -6.0, -41.0, -9.0, -3.0, 0.0, 0.0, -6.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -96.0, 33.0, -12.0, 228.0, -23.0, -21.0, 0.0, 0.0, -12.0, -2.0, -4.0, 4.0, -1.0, 0.0, 1.0, 0.0, -329.0, -22.0, -34.0, -726.0, -147.0, -21.0, 0.0, 0.0, -2.0, 4.0, -1.0, 0.0, 2.0, -7.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.0, 88.0, -162.0, -19.0, -11.0, 21.0, 31.0, 37.0, -31.0, 53.0, -5.0, -15.0, -3.0, -11.0, 9.0, 3.0, 0.0, 0.0, -2.0, 0.0, 1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0, -162.0, -102.0, -37.0, 30.0, 19.0, 23.0, -18.0, 9.0, 1.0, -6.0, -6.0, 22.0, -2.0, 3.0, 1.0, -2.0, 0.0, -1.0, 26.0, -25.0, 66.0, 52.0, -641.0, -153.0, -13.0, -9.0, 2.0, -3.0, -29.0, 8.0, -6.0, -2.0, 0.0, -6.0, 2.0, -4.0, 1.0, 0.0, -26.0, -11.0, -1.0, -10.0, -6.0, -13.0, 66.0, -1337.0, -879.0, -207.0, 1.0, -1.0, 8.0, -30.0, -24.0, -18.0, -16.0, 1.0, 9.0, 1.0, -24.0, -8.0, 9.0, -17.0, -13.0, 75.0, 19.0, -8.0, -29.0, 24.0, 0.0, 0.0, -1.0, 1.0, -25.0, 36.0, -7.0, -22.0, 0.0, -3.0, 1.0, -1.0, 187.0, -46.0, -6.0, 74.0, 5.0, -10.0, -5.0, -4.0, -16.0, 10.0, -5.0, -5.0, 2.0, -4.0, 5.0, -2.0, -2.0, 1.0, -1.0, 0.0, -16.0, -12.0, 1.0, -13.0, -17.0, -111.0, -186.0, 73.0, -1.0, -2.0, -277.0, -77.0, -27.0, 106.0, 16.0, 5.0, -12.0, -15.0, -13.0, -30.0, -1.0, 1.0, 0.0, 36.0, -10.0, 4.0, 607.0, 262.0, 533.0, -1530.0, -1630.0, 304.0, 8.0, -6.0, 1.0, 1.0, 0.0, -1.0, 5.0, -2.0, 0.0, -1.0, -1.0, -4.0, -44.0, -22.0, -64.0, -46.0, 537.0, 430.0, 268.0, -1553.0, -2040.0, -486.0, -3.0, -23.0, 20.0, 41.0, -1.0, 2.0, -21.0, -4.0, -1.0, -3.0, -84.0, 50.0, -177.0, 26.0, 5.0, -12.0, 2.0, -4.0, 7.0, 1.0, -115.0, -305.0, -310.0, 138.0, -186.0, 246.0, -96.0, 17.0, 0.0, 0.0, 4.0, -2.0, 1.0, 1.0, -3.0, 2.0, -1.0, 0.0, -15.0, 68.0, 0.0, 2.0, -3.0, 0.0, -5.0, 0.0, -1.0, 1.0, -5.0, 6.0, 0.0, 0.0, 0.0, 0.0, -235.0, -98.0, -2.0, 2.0, 9.0, -40.0, -1.0, -2.0, -33.0, -9.0, -5.0, -4.0, 5662.0, -3849.0, 1941.0, -124.0, 210.0, 160.0, -24721.0, -72945.0, 4099.0, -21914.0, 1345.0, -555.0, 23637393.0, -5516830.0, 17737677.0, 43330654.0, -44668315.0, 14540723.0, -824.0, -2086.0, -4423.0, -41801.0, 5562.0, -11664.0, 960.0, -125.0, 2001.0, -149.0, 587.0, -350.0, 23.0, -52.0, -3.0, 3.0, -248.0, -148.0, -40.0, 86.0, 2.0, 0.0, 21.0, -82.0, 11.0, 8.0, -8.0, 0.0, -30.0, -33.0, -22.0, 46.0, 0.0, -191.0, -168.0, -135.0, 27.0, -85.0, 14.0, 232.0, 217.0, 59.0, 5.0, 12.0, -5.0, 2.0, -24.0, -26.0, -52.0, 13.0, -3.0, 18.0, 26.0, 45.0, 32.0, -169.0, 14.0, -6.0, -3.0, 4.0, -5.0, 2.0, 6.0, 2.0, -2.0, 3.0, 20.0, -15.0, 0.0, 10.0, -486.0, -8.0, 4.0, -114.0, 102.0, -188.0, 23.0, -67.0, 6.0, 12.0, -43.0, -1.0, -32.0, 2.0, 15.0, 9.0, 16.0, -36.0, -6.0, -2.0, 14.0, -5.0, 17.0, -15.0, -28.0, 307.0, 289.0, 69.0, 2.0, -7.0, 3.0, -1.0, -1.0, 1.0, -16.0, -811.0, 287.0, -68.0, 0.0, 0.0, 0.0, -1.0, 16.0, -7.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, -3.0, -4.0, 2.0, 3.0, -29.0, 34.0, 59.0, -15.0, -3.0, -3.0, -1.0, 0.0, -2.0, -3.0, 3.0, -19.0, 0.0, 0.0, 0.0, 0.0, -15.0, 1.0, 5.0, 2.0, 0.0, 0.0, -1.0, -5.0, 0.0, -1.0, -120.0, 84.0, 7.0, -30.0, -7.0, -3.0, -1.0, 0.0, 0.0, -1.0, 9.0, -6.0, -186.0, -11.0, 13.0, -57.0, 1.0, 4.0, 1.0, -1.0, 0.0, 0.0, -5.0, 796.0, 46.0, 5.0, -1.0, -6.0, -10.0, 228.0, 5.0, -6.0, 1.0, -5.0, 0.0, 0.0, -6.0, -2.0, 148.0, 137.0, 10.0, 28.0, 430546.0, -279834.0, 488902.0, 664558.0, -746515.0, 243112.0, -39.0, -37.0, -13.0, -174.0, 6.0, -25.0, 2.0, -3.0, -4.0, -2.0, 0.0, 4.0, -5.0, 70.0, 82.0, 20.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, -27.0, 430.0, 226.0, -53.0, 1.0, 1.0, 0.0, 1.0, 1.0, -7.0, 2.0, 1.0, -3.0, -8.0, 1.0, 0.0, -1.0, 12.0, -2.0, -5.0, 4.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 9.0, 33.0, 4.0, 0.0, 0.0, 0.0, -321.0, 4.0, 1.0, 0.0, 0.0, 1.0, 0.0, 106.0, -22.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 7006.0, -9443.0, 12833.0, 11137.0, -14037.0, 4575.0, -2.0, 0.0, -1.0, -6.0, 1.0, 1.0, 4.0, 6.0, 16.0, 2.0, 55.0, -10.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, -4.0, -2.0, 0.0, -351.0, 24.0, 0.0, 0.0, 8.0, 1.0, 30.0, -5.0, -12.0, 10.0, -4.0, 1.0, -1.0, -2.0, 0.0, 0.0, 4.0, 0.0, 17.0, -3.0, 0.0, -2.0, 2.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0) }
    override val tabr by lazy { doubleArrayOf(14575.0, -26192.0, -144864.0, 2.0, -22.0, 15.0, -8.0, -21.0, -148.0, -104.0, -14.0, -75.0, 15.0, 2.0, -5.0, -3.0, -1.0, 0.0, 0.0, 0.0, 0.0, 21.0, -2.0, 7.0, -5.0, -3.0, 0.0, 0.0, 83.0, -94.0, 9.0, -67.0, -29.0, 50.0, 1.0, -1.0, 3.0, 2.0, 0.0, 0.0, 4.0, 3.0, 1.0, 1.0, -1.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, 1.0, 0.0, 0.0, -2.0, 3.0, -5.0, -2.0, -1.0, 1.0, 197.0, 511.0, -82.0, 189.0, -28.0, -12.0, 0.0, 0.0, 0.0, -1.0, 6.0, -1.0, 0.0, 1.0, 30.0, -30.0, -37.0, -25.0, 6.0, 21.0, 0.0, 0.0, 16.0, -139.0, 43.0, -28.0, 4.0, 6.0, 0.0, 3.0, 4.0, 0.0, 1.0, 1.0, -13.0, 0.0, -4.0, 0.0, 0.0, 1.0, 150.0, 135.0, -291.0, 436.0, -560.0, -343.0, 1.0, 3.0, 8.0, -15.0, -13.0, -5.0, -17.0, -3.0, 1.0, -6.0, -314.0, 428.0, 606.0, 758.0, 1230.0, -411.0, 0.0, -1.0, 11.0, -14.0, 4.0, 1.0, 221.0, 157.0, 1.0, 132.0, -25.0, 3.0, 12.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1487.0, -108.0, 707.0, -79.0, -950.0, -190.0, 177.0, 582.0, -676.0, 399.0, -281.0, -396.0, 0.0, 52.0, 39.0, -130.0, 2.0, 1.0, 12.0, 148.0, -34.0, 23.0, 1.0, 27.0, -20.0, 13.0, 1.0, -1.0, 198.0, -34.0, -21.0, -80.0, -99.0, 332.0, -307.0, 9.0, -15.0, -125.0, 330.0, -231.0, 236.0, 139.0, -36.0, 74.0, 7.0, 3.0, -588.0, -1722.0, 3623.0, -1245.0, 187.0, 4366.0, -72.0, -75.0, 11.0, -33.0, 174.0, -467.0, 444.0, 9.0, 11.0, 180.0, -6.0, -39.0, 8.0, -7.0, -126.0, -500.0, 599.0, -317.0, 224.0, 355.0, -590.0, -39.0, 134.0, -379.0, -7.0, -3.0, 494.0, -628.0, 893.0, 490.0, 712.0, -7.0, -7.0, 35.0, -720.0, 50.0, -321.0, 72.0, 443.0, 106.0, 74.0, 82.0, 112.0, -84.0, -6.0, -3.0, 5.0, 4.0, 58.0, 7.0, -2.0, 38.0, 6.0, 92.0, -20.0, 14.0, 33.0, 13.0, -11189.0, -2.0, -11.0, -8.0, 106.0, -35.0, 58.0, 52.0, 132.0, 170.0, -32.0, 63.0, -2.0, -6.0, 6.0, 7.0, -1.0, 6.0, 452.0, 155.0, 9.0, 209.0, 788.0, -318.0, 511.0, 616.0, -5.0, 3.0, 142.0, 303.0, -280.0, 32.0, 21.0, -69.0, 984.0, -291.0, 340.0, 562.0, 30.0, 2.0, 171.0, -51.0, 27.0, 28.0, -1570.0, -2053.0, 3702.0, -3593.0, 4012.0, 5467.0, -14.0, -9.0, -1.0, -6.0, 4.0, -11.0, 0.0, -9.0, -11.0, 0.0, 15.0, 2.0, 1133.0, -1366.0, 1961.0, 1134.0, -867.0, -3010.0, 6041.0, -2049.0, 142.0, 7138.0, -79.0, -103.0, 73.0, -18.0, 17.0, -9.0, 79.0, -372.0, 97.0, -17.0, 182.0, -118.0, 33577.0, -675.0, -99.0, -347.0, -91.0, -11.0, 1.0, -7.0, 158.0, -111.0, 136.0, 474.0, 50.0, 16.0, 9739.0, 51.0, 19.0, 2.0, -105.0, 36.0, -20.0, -11.0, -12.0, 56.0, -13.0, 2.0, 2030.0, -22.0, 11.0, -4.0, 9.0, -6.0, 1.0, 33.0, 10.0, 18.0, 4.0, 21.0, 53.0, 3.0, 19.0, 2.0, 130.0, 37.0, -147.0, -54.0, -22.0, 6.0, 7.0, -19.0, 22.0, 36.0, -4.0, 8.0, -949.0, 2911.0, -2221.0, -697.0, 371.0, -808.0, 1203.0, -1117.0, 191.0, 189.0, -549.0, 77.0, -321201.0, 19788.0, 2622593.0, 5990.0, 8667033.0, 114.0, -34455835.0, 86.0, -92.0, -493.0, 179.0, 807.0, -499.0, 281.0, 225.0, 51.0, -34.0, -88.0, -30.0, 13.0, 0.0, -3.0, 20.0, -43.0, 201.0, 33.0, -93.0, 2.0, -1034.0, -518.0, -63.0, -545.0, 104.0, -20.0, 43.0, -64.0, -34.0, -7.0, 0.0, 0.0, -61.0, -159.0, -143.0, -8.0, -392.0, -157.0, -204.0, 309.0, -24.0, -248.0, 55.0, -40.0, -6.0, 91.0, -16.0, 57.0, -41.0, 18.0, 197.0, -20.0, -668.0, -150.0, -192.0, -216.0, 39.0, -84.0, -62.0, -59.0, -4.0, 8.0, -7.0, -1.0, -352.0, -100.0, -10.0, -158.0, 61.0, 55.0, 32493.0, -49.0, 107.0, 344.0, -395.0, 227.0, -154.0, -238.0, 123.0, 104.0, -205.0, 348.0, -449.0, -236.0, -54.0, -19.0, -6.0, 21.0, -790.0, 27.0, -5.0, 30.0, -846.0, 154.0, -26.0, -920.0, 0.0, 12.0, -54.0, 21.0, 11.0, -10.0, 137.0, 132.0, 109.0, -337.0, -11.0, 2.0, -272.0, 467.0, -511.0, 179.0, -1.0, -8.0, 7.0, -32.0, -44.0, 2.0, -26.0, 101.0, -32.0, 6.0, -98.0, 48.0, -42.0, -53.0, -1222.0, 1601.0, -3775.0, 656.0, 83.0, 46.0, 16.0, -28.0, 0.0, 7.0, -66.0, -14.0, -6.0, 13.0, 12.0, 0.0, -58.0, 91.0, -123.0, -58.0, -12.0, 4.0, -114.0, 423.0, -111.0, 12.0, 112.0, 27.0, -19072.0, 71.0, 100.0, 410.0, 107.0, 15.0, 24.0, 3.0, -214.0, 30.0, 49.0, 44.0, 5017.0, -27.0, 167.0, -80.0, 8.0, 153.0, 4.0, 7.0, -219.0, -35.0, 244.0, 694.0, -762.0, 2.0, -84.0, -49.0, -28.0, 158.0, -4.0, 56.0, -14.0, 0.0, 9.0, 12.0, 7.0, 18.0, 2.0, -7.0, -15426.0, 91.0, 25800.0, -15.0, 144767.0, -53.0, -287824.0, -24.0, 19.0, -9.0, 6.0, 7.0, 0.0, 0.0, -3.0, 8.0, -5.0, -3.0, -232.0, 53.0, -1.0, -271.0, 4.0, -12.0, -8.0, 30.0, -8.0, -25.0, -253.0, -150.0, -105.0, 470.0, -37.0, -29.0, -59.0, -6.0, -24.0, -5.0, 9.0, -18.0, 1784.0, 3.0, -54.0, 13.0, -12.0, 7.0, -116.0, 144.0, -353.0, 52.0, -4.0, -12.0, -17.0, -14.0, -1340.0, -64.0, 10.0, -116.0, -24.0, -2.0, 190.0, 131.0, 130.0, -307.0, -1.0, 9.0, 5.0, -7.0, -10.0, 56.0, -33.0, 0.0, -14.0, 3.0, 2.0, -12.0, -635.0, -160.0, 64.0, -44.0, 2712.0, -3.0, -3606.0, -1.0, 774.0, 1.0, 133.0, -1.0, -19.0, 0.0, 5.0, -59.0, -5.0, 14.0, -45.0, 5.0, -140.0, -8.0, 15.0, -28.0, 379.0, 1.0, 6.0, 3.0, 55.0, 0.0, -54.0, 0.0, 3.0, -33.0, -3.0, 4.0, 0.0, -4.0, 0.0, -1.0, 200.0, 0.0, -17.0, -1.0, 2.0, -20.0, -2.0, 0.0, 111.0, 0.0, 1.0, -12.0, 64.0, 0.0, 38.0, 0.0, 23.0, 0.0, 3.0, 13.0) }
    override val args by lazy { intArrayOf(0, 3, 3, 4, 3, -8, 4, 3, 5, 1, 2, 2, 5, -5, 6, 2, 4, 4, 3, -8, 4, 5, 5, -5, 6, 1, 3, 2, 2, 1, 3, -8, 4, 0, 3, 3, 2, -7, 3, 4, 4, 2, 3, 7, 3, -13, 4, -1, 5, 0, 2, 8, 2, -13, 3, 2, 3, 1, 3, -2, 4, 2, 6, 0, 3, 1, 2, -8, 3, 12, 4, 1, 3, 6, 2, -10, 3, 3, 5, 1, 1, 1, 7, 0, 2, 1, 5, -2, 6, 1, 2, 1, 5, -3, 6, 0, 3, 1, 3, -2, 4, 1, 5, 0, 3, 3, 3, -6, 4, 2, 5, 1, 3, 1, 1, -5, 2, 4, 3, 0, 2, 8, 3, -15, 4, 2, 3, 4, 3, -7, 4, -3, 5, 0, 3, 2, 2, -7, 3, 7, 4, 0, 2, 2, 5, -4, 6, 1, 1, 1, 6, 2, 2, 2, 5, -6, 6, 0, 2, 9, 3, -17, 4, 2, 3, 3, 2, -5, 3, 1, 5, 0, 3, 2, 3, -4, 4, 2, 5, 0, 3, 2, 3, -4, 4, 1, 5, 0, 3, 3, 2, -5, 3, 2, 5, 0, 2, 1, 5, -1, 6, 0, 3, 3, 2, -6, 3, 2, 4, 0, 2, 1, 3, -2, 4, 2, 2, 2, 5, -3, 6, 0, 1, 2, 6, 1, 2, 3, 5, -5, 6, 1, 1, 1, 5, 2, 3, 4, 3, -8, 4, 2, 5, 0, 2, 1, 5, -5, 6, 1, 2, 7, 3, -13, 4, 2, 2, 2, 5, -2, 6, 0, 2, 10, 3, -19, 4, 0, 2, 3, 5, -4, 6, 0, 2, 3, 2, -5, 3, 2, 2, 2, 3, -4, 4, 2, 2, 5, 2, -8, 3, 1, 2, 3, 5, -3, 6, 0, 2, 6, 3, -11, 4, 1, 2, 1, 1, -4, 3, 1, 2, 4, 5, -5, 6, 0, 1, 2, 5, 1, 2, 3, 3, -6, 4, 2, 2, 5, 3, -9, 4, 2, 2, 6, 2, -10, 3, 0, 2, 2, 2, -3, 3, 2, 2, 4, 3, -8, 4, 1, 2, 4, 3, -7, 4, 2, 2, 5, 3, -10, 4, 1, 2, 3, 3, -5, 4, 2, 2, 1, 2, -2, 3, 1, 2, 7, 2, -11, 3, 0, 2, 2, 3, -3, 4, 1, 2, 1, 3, -1, 4, 0, 2, 4, 2, -7, 3, 0, 2, 4, 2, -6, 3, 2, 1, 1, 4, 1, 2, 8, 3, -14, 4, 0, 2, 1, 3, -5, 5, 0, 2, 1, 3, -3, 4, 1, 2, 7, 3, -12, 4, 1, 2, 1, 2, -1, 3, 1, 2, 2, 3, -5, 4, 0, 2, 1, 3, -4, 5, 1, 2, 6, 3, -10, 4, 1, 2, 3, 3, -7, 4, 0, 3, 1, 3, -4, 5, 2, 6, 0, 3, 1, 3, -1, 5, -5, 6, 0, 2, 5, 3, -8, 4, 1, 2, 1, 3, -3, 5, 1, 3, 1, 3, -5, 5, 5, 6, 0, 2, 2, 2, -4, 3, 1, 2, 6, 2, -9, 3, 0, 2, 4, 3, -6, 4, 1, 3, 1, 3, -3, 5, 2, 6, 0, 2, 1, 3, -5, 6, 1, 2, 1, 3, -2, 5, 2, 3, 1, 3, -4, 5, 5, 6, 0, 3, 1, 3, -1, 5, -2, 6, 0, 3, 1, 3, -3, 5, 3, 6, 0, 2, 1, 3, -4, 6, 0, 3, 1, 3, -2, 5, 1, 6, 0, 2, 5, 2, -9, 3, 0, 2, 3, 3, -4, 4, 1, 2, 3, 2, -4, 3, 2, 2, 1, 3, -3, 6, 1, 3, 1, 3, -2, 5, 2, 6, 0, 3, 1, 3, 1, 5, -5, 6, 1, 2, 1, 3, -1, 5, 1, 3, 1, 3, -3, 5, 5, 6, 1, 3, 1, 3, 2, 5, -7, 6, 0, 2, 1, 3, -2, 6, 1, 2, 2, 3, -2, 4, 1, 3, 3, 2, -4, 3, 1, 5, 0, 2, 10, 3, -17, 4, 1, 3, 1, 3, 2, 5, -6, 6, 1, 2, 1, 3, -1, 6, 0, 3, 1, 3, -2, 5, 4, 6, 0, 2, 7, 3, -15, 4, 0, 2, 1, 3, -2, 7, 0, 3, 1, 3, 1, 5, -3, 6, 0, 2, 1, 3, -2, 8, 0, 2, 1, 3, -1, 7, 0, 2, 1, 3, -1, 8, 0, 2, 8, 2, -14, 3, 1, 3, 3, 2, -8, 3, 4, 4, 1, 3, 1, 3, 4, 5, -10, 6, 1, 3, 1, 3, 2, 5, -5, 6, 2, 3, 5, 3, -8, 4, 3, 5, 2, 1, 1, 12, 3, 3, 3, 3, -8, 4, 3, 5, 2, 3, 1, 3, -2, 5, 5, 6, 2, 3, 3, 2, -6, 3, 4, 4, 0, 2, 8, 2, -12, 3, 1, 3, 1, 3, 1, 5, -2, 6, 0, 2, 9, 3, -15, 4, 2, 2, 1, 3, 1, 6, 1, 2, 1, 10, -1, 11, 0, 1, 2, 4, 1, 2, 1, 3, 1, 5, 1, 2, 8, 3, -13, 4, 1, 2, 3, 2, -6, 3, 0, 2, 1, 3, -4, 4, 1, 2, 5, 2, -7, 3, 1, 2, 7, 3, -11, 4, 1, 2, 1, 1, -3, 3, 0, 2, 1, 3, 2, 5, 0, 2, 2, 3, -6, 4, 0, 2, 6, 3, -9, 4, 1, 2, 2, 2, -2, 3, 1, 2, 5, 3, -7, 4, 2, 2, 4, 3, -5, 4, 2, 2, 1, 2, -3, 3, 0, 2, 7, 2, -10, 3, 0, 2, 3, 3, -3, 4, 0, 2, 2, 3, -1, 4, 0, 2, 4, 2, -5, 3, 1, 2, 1, 3, 1, 4, 0, 2, 2, 3, -5, 5, 0, 2, 8, 3, -12, 4, 0, 1, 1, 2, 1, 3, 2, 3, -5, 5, 2, 6, 0, 2, 2, 3, -4, 5, 1, 3, 2, 3, -6, 5, 5, 6, 0, 2, 7, 3, -10, 4, 0, 3, 2, 3, -4, 5, 2, 6, 0, 3, 2, 3, -1, 5, -5, 6, 1, 2, 6, 3, -8, 4, 1, 2, 2, 3, -3, 5, 1, 3, 2, 3, -5, 5, 5, 6, 1, 2, 2, 2, -5, 3, 0, 2, 6, 2, -8, 3, 0, 3, 2, 3, -4, 5, 3, 6, 0, 3, 2, 3, -3, 5, 1, 6, 0, 2, 5, 3, -6, 4, 1, 3, 2, 3, -3, 5, 2, 6, 0, 2, 2, 3, -5, 6, 1, 2, 2, 3, -2, 5, 1, 3, 2, 3, -4, 5, 5, 6, 1, 2, 2, 3, -4, 6, 0, 2, 4, 3, -4, 4, 0, 2, 3, 2, -3, 3, 1, 2, 2, 3, -3, 6, 1, 3, 2, 3, -2, 5, 2, 6, 0, 2, 2, 3, -1, 5, 1, 2, 2, 3, -2, 6, 0, 2, 3, 3, -2, 4, 1, 2, 2, 3, -1, 6, 0, 2, 2, 3, -2, 7, 0, 3, 2, 3, 2, 5, -5, 6, 0, 3, 6, 3, -8, 4, 3, 5, 1, 1, 2, 12, 3, 3, 2, 3, -8, 4, 3, 5, 1, 3, 2, 3, -2, 5, 5, 6, 0, 2, 8, 2, -11, 3, 0, 2, 2, 3, 1, 5, 0, 2, 5, 2, -6, 3, 1, 2, 8, 3, -11, 4, 0, 2, 1, 1, -2, 3, 0, 2, 7, 3, -9, 4, 0, 2, 2, 2, -1, 3, 1, 2, 6, 3, -7, 4, 0, 2, 5, 3, -5, 4, 0, 2, 7, 2, -9, 3, 0, 2, 4, 3, -3, 4, 0, 2, 4, 2, -4, 3, 0, 2, 3, 3, -5, 5, 0, 2, 1, 2, 1, 3, 0, 2, 3, 3, -4, 5, 1, 2, 8, 3, -10, 4, 0, 2, 7, 3, -8, 4, 0, 2, 3, 3, -3, 5, 0, 2, 6, 2, -7, 3, 0, 2, 6, 3, -6, 4, 0, 2, 3, 3, -2, 5, 1, 2, 3, 3, -4, 6, 0, 2, 5, 3, -4, 4, 0, 2, 3, 2, -2, 3, 0, 2, 3, 3, -3, 6, 0, 2, 3, 3, -1, 5, 0, 2, 3, 3, -2, 6, 0, 1, 3, 12, 3, 2, 5, 2, -5, 3, 0, 2, 1, 1, -1, 3, 0, 1, 2, 2, 0, 2, 7, 2, -8, 3, 0, 2, 4, 2, -3, 3, 0, 2, 4, 3, -5, 5, 0, 2, 4, 3, -4, 5, 0, 2, 4, 3, -3, 5, 0, 2, 6, 2, -6, 3, 0, 2, 4, 3, -2, 5, 0, 1, 4, 12, 1, 2, 8, 2, -9, 3, 0, 2, 5, 2, -4, 3, 0, 1, 1, 1, 0, 2, 7, 2, -7, 3, 1, 2, 5, 3, -5, 5, 0, 2, 9, 2, -10, 3, 0, 2, 6, 2, -5, 3, 0, 2, 8, 2, -8, 3, 0, 2, 10, 2, -11, 3, 0, 2, 9, 2, -9, 3, 0, 2, 10, 2, -10, 3, 0, 2, 11, 2, -11, 3, 0, 2, 2, 1, -1, 3, 0, -1) }

    override val maxargs = 18
    override val max_harmonic by lazy { intArrayOf(2, 11, 14, 19, 6, 10, 2, 2, 0, 1, 1, 4, 0, 0, 0, 0, 0, 0) }
    override val max_power_of_t = 3
    override val distance = 1.000139872959708


}