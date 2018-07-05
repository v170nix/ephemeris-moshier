package net.arwix.astronomy2.ephemeris.moshier

import net.arwix.astronomy2.core.ARCSEC_TO_RAD
import net.arwix.astronomy2.core.JT
import net.arwix.astronomy2.core.ephemeris.coordinates.HeliocentricEclipticEphemeris
import net.arwix.astronomy2.core.ephemeris.precession.ID_PRECESSION_DE4xx
import net.arwix.astronomy2.core.ephemeris.precession.ID_PRECESSION_WILLIAMS_1994
import net.arwix.astronomy2.core.ephemeris.precession.createEclipticPrecessionMatrix
import net.arwix.astronomy2.core.kepler.ID_EARTH_KEPLER_ELEMENTS
import net.arwix.astronomy2.core.kepler.getAscendingNode
import net.arwix.astronomy2.core.kepler.getLongitudeOfMoon
import net.arwix.astronomy2.core.kepler.getSimonJ2000KeplerElements
import net.arwix.astronomy2.core.vector.Matrix
import net.arwix.astronomy2.core.vector.Matrix.Companion.AXIS_X
import net.arwix.astronomy2.core.vector.Matrix.Companion.AXIS_Z
import net.arwix.astronomy2.core.vector.RectangularVector
import net.arwix.astronomy2.core.vector.Vector
import net.arwix.astronomy2.core.vector.convert
import kotlin.math.abs
import kotlin.math.atan2
import kotlin.math.cos
import kotlin.math.sin


fun getMoshierEphemeris(id: IdMoshierBody): HeliocentricEclipticEphemeris = when (id) {
    ID_MOSHIER_SUN -> SunMoshierEphemeris
    ID_MOSHIER_MERCURY -> MercuryMoshierEphemeris
    ID_MOSHIER_VENUS -> VenusMoshierEphemeris
    ID_MOSHIER_EARTH -> EarthMoshierEphemeris
    ID_MOSHIER_LIBRATION -> LibrationMoshierEphemeris
    ID_MOSHIER_EM_BARYCENTER -> EarthMoonBarycenterMoshierEphemeris
    ID_MOSHIER_MOON -> EarthMoonBarycenterMoshierEphemeris
    ID_MOSHIER_MARS -> MarsMoshierEphemeris
    ID_MOSHIER_JUPITER -> JupiterMoshierEphemeris
    ID_MOSHIER_SATURN -> SaturnMoshierEphemeris
    ID_MOSHIER_URANUS -> UranusMoshierEphemeris
    ID_MOSHIER_NEPTUNE -> NeptuneMoshierEphemeris
    ID_MOSHIER_PLUTO -> PlutoMoshierEphemeris
    else -> throw IndexOutOfBoundsException()
}




object SunMoshierEphemeris : HeliocentricEclipticEphemeris {
    override fun getCoordinates(jT: JT) = RectangularVector()
}
object MercuryMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(MercuryMoshierData)
object VenusMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(VenusMoshierData)
object MarsMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(MarsMoshierData)
object JupiterMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(JupiterMoshierData)
object SaturnMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(SaturnMoshierData)
object UranusMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(UranusMoshierData)
object NeptuneMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(NeptuneMoshierData)
object PlutoMoshierEphemeris: HeliocentricEclipticEphemeris by BaseImpl(PlutoMoshierData)

object EarthMoonBarycenterMoshierEphemeris: HeliocentricEclipticEphemeris {
    override fun getCoordinates(jT: JT): Vector = RectangularVector(
            g3plan(jT, EarthMoonBarycenterMoshierData.args, EarthMoonBarycenterMoshierData.distance,
                    EarthMoonBarycenterMoshierData.tabb, EarthMoonBarycenterMoshierData.tabl,
                    EarthMoonBarycenterMoshierData.tabr, EarthMoonBarycenterMoshierData.max_harmonic,
                    EarthMoonBarycenterMoshierData.max_power_of_t, EarthMoonBarycenterMoshierData.maxargs,
                    EarthMoonBarycenterMoshierData.timescale, EarthMoonBarycenterMoshierData.trunclvl, false))

}

object MoonMoshierEphemeris: HeliocentricEclipticEphemeris {
    override fun getCoordinates(jT: JT): Vector {
        val moonLat = g1plan(jT, MoonLatMoshierData.args, MoonLatMoshierData.tabl, MoonLatMoshierData.max_harmonic, MoonLatMoshierData.timescale)
        val moon = g2plan(jT, MoonLonMoshierData.args, MoonLonMoshierData.distance,
                MoonLonMoshierData.tabl, MoonLonMoshierData.tabr, MoonLonMoshierData.max_harmonic, MoonLonMoshierData.timescale, moonLat)
        val precession = createEclipticPrecessionMatrix(ID_PRECESSION_WILLIAMS_1994, jT).transpose()
        return precession * moon
    }
}

object EarthMoshierEphemeris: HeliocentricEclipticEphemeris {
    override fun getCoordinates(jT: JT): Vector {
        val p = EarthMoonBarycenterMoshierEphemeris.getCoordinates(jT)
        val moon = MoonMoshierEphemeris.getCoordinates(jT)
        val earthMoonRatio = 2.7068700387534E7 / 332946.050895
        return p - moon * (1.0 / (earthMoonRatio + 1.0))
    }
}

object LibrationMoshierEphemeris: HeliocentricEclipticEphemeris {
    override fun getCoordinates(jT: JT): Vector {
        val p = convert<RectangularVector>(g3plan(jT, LibrationMoshierData.args, LibrationMoshierData.distance,
                LibrationMoshierData.tabb, LibrationMoshierData.tabl,
                LibrationMoshierData.tabr, LibrationMoshierData.max_harmonic,
                LibrationMoshierData.max_power_of_t, LibrationMoshierData.maxargs,
                LibrationMoshierData.timescale, LibrationMoshierData.trunclvl, true))
        p.x = -getSimonJ2000KeplerElements(ID_EARTH_KEPLER_ELEMENTS).getLongitude(jT) - 0.047 * ARCSEC_TO_RAD

        val LP_equinox = getLongitudeOfMoon(jT)
        val NF_arcsec = getAscendingNode(jT)

        // phi+psi
        p[2] += LP_equinox + 6.48e5 * ARCSEC_TO_RAD
        if (p[2] < -6.45e5 * ARCSEC_TO_RAD) p[2] += 1.296e6 * ARCSEC_TO_RAD
        if (p[2] > 6.45e5 * ARCSEC_TO_RAD) p[2] -= 1.296e6 * ARCSEC_TO_RAD

        // phi
        p[0] += LP_equinox - NF_arcsec + 6.48e5 * ARCSEC_TO_RAD
        if (p[0] < -6.45e5 * ARCSEC_TO_RAD) p[0] += 1.296e6 * ARCSEC_TO_RAD
        if (p[0] > 6.45e5 * ARCSEC_TO_RAD) p[0] -= 1.296e6 * ARCSEC_TO_RAD
        p[2] -= p[0]

        // From Euler angles to matrix M
        val x = Matrix(AXIS_Z, p[2])
        val y = Matrix(AXIS_X, p[1])
        val z = Matrix(AXIS_Z, p[0])
        val mM = x * y * z

        // Get rotation matrix around x axis with eps
        val mQ = Matrix(AXIS_X, 84381.406173 * ARCSEC_TO_RAD)

        // Get precession matrix
        var mP = createEclipticPrecessionMatrix(ID_PRECESSION_DE4xx, jT)

        // Precess Q
        mP = mP * mQ

        // Space to body
        val mM2000 = mM * mP

        // Get back the Euler angles, now equatorial (as JPL ones)
        val M = mM2000.toArray() // mM2000.array
        val phi = atan2(M[2][0], -M[2][1])
        var a = M[0][2]
        val b = M[1][2]
        /* psi = zatan2( b, a ); */
        val psi = atan2(a, b)

        if (abs(a) > abs(b)) a /= sin(psi) else a = b / cos(psi)

        /* theta = zatan2( M[2][2], a ); */
        val theta = atan2(a, M[2][2])

        p[0] = phi
        p[1] = theta
        p[2] = psi
        return RectangularVector(p)
    }
}

private class BaseImpl(val data: MoshierData) : HeliocentricEclipticEphemeris {

    override fun getCoordinates(jT: JT): Vector = RectangularVector(
            gplan(jT, data.args, data.distance, data.tabb, data.tabl,
                    data.tabr, data.max_harmonic, data.max_power_of_t,
                    data.maxargs, data.timescale, data.trunclvl))
}
