package net.arwix.astronomy2.ephemeris.moshier

import kotlinx.coroutines.experimental.CommonPool
import kotlinx.coroutines.experimental.async
import net.arwix.astronomy2.core.*
import net.arwix.astronomy2.core.ephemeris.coordinates.getCoroutineHeliocentricEclipticCoordinates
import net.arwix.astronomy2.core.ephemeris.coordinates.getHeliocentricEclipticCoordinates
import net.arwix.astronomy2.core.ephemeris.precession.*
import net.arwix.astronomy2.core.kepler.ID_EARTH_KEPLER_ELEMENTS
import net.arwix.astronomy2.core.kepler.getAscendingNode
import net.arwix.astronomy2.core.kepler.getLongitudeOfMoon
import net.arwix.astronomy2.core.kepler.getSimonJ2000KeplerElements
import net.arwix.astronomy2.core.vector.Matrix
import net.arwix.astronomy2.core.vector.RectangularVector
import net.arwix.astronomy2.core.vector.Vector
import net.arwix.astronomy2.core.vector.convert
import kotlin.math.abs
import kotlin.math.atan2
import kotlin.math.cos
import kotlin.math.sin


@Heliocentric @Ecliptic @J2000
fun createMoshierCoordinates(idMoshierBody: IdMoshierBody): getHeliocentricEclipticCoordinates {
    val data = when (idMoshierBody) {
        ID_MOSHIER_SUN -> return {_
            ->RectangularVector()
        }
        ID_MOSHIER_MERCURY -> MercuryMoshierData
        ID_MOSHIER_VENUS -> VenusMoshierData
        ID_MOSHIER_EM_BARYCENTER -> return {jT ->
            g3plan(jT, EarthMoonBarycenterMoshierData.args, EarthMoonBarycenterMoshierData.distance,
                    EarthMoonBarycenterMoshierData.tabb, EarthMoonBarycenterMoshierData.tabl,
                    EarthMoonBarycenterMoshierData.tabr, EarthMoonBarycenterMoshierData.max_harmonic,
                    EarthMoonBarycenterMoshierData.max_power_of_t, EarthMoonBarycenterMoshierData.maxargs,
                    EarthMoonBarycenterMoshierData.timescale, EarthMoonBarycenterMoshierData.trunclvl, false)

        }
        ID_MOSHIER_MOON -> return {jT -> createMoonMoshierCoordinates(jT).invoke(jT) }
        ID_MOSHIER_EARTH -> return {jT ->
            createEarthMoshierCoordinates(
                    createMoshierCoordinates(ID_MOSHIER_EM_BARYCENTER)
            ) { _ -> createMoonMoshierCoordinates(jT).invoke(jT) }.invoke(jT)
        }
        ID_MOSHIER_LIBRATION -> return createLibrationMoshierCoordinates()
        ID_MOSHIER_MARS -> MarsMoshierData
        ID_MOSHIER_JUPITER -> JupiterMoshierData
        ID_MOSHIER_SATURN -> SaturnMoshierData
        ID_MOSHIER_URANUS -> UranusMoshierData
        ID_MOSHIER_NEPTUNE -> NeptuneMoshierData
        ID_MOSHIER_PLUTO -> PlutoMoshierData
        else -> throw IndexOutOfBoundsException()
    }
    return { jT -> getMoshierVector(jT, data) }
}

@Heliocentric @Ecliptic @J2000
fun createMoonMoshierCoordinates(jT0: JT,
                                 @Ecliptic precessionElements: PrecessionElements =
                                         createPrecessionElements(ID_PRECESSION_WILLIAMS_1994, jT0)
): getHeliocentricEclipticCoordinates {
    return {jT ->
        val moonLat = g1plan(jT, MoonLatMoshierData.args, MoonLatMoshierData.tabl, MoonLatMoshierData.max_harmonic, MoonLatMoshierData.timescale)
        val moon = g2plan(jT, MoonLonMoshierData.args, MoonLonMoshierData.distance,
                MoonLonMoshierData.tabl, MoonLonMoshierData.tabr, MoonLonMoshierData.max_harmonic, MoonLonMoshierData.timescale, moonLat)
        precessionElements.transformToJ2000(moon)
    }
}

@Heliocentric @Ecliptic @J2000
inline fun createEarthMoshierCoordinates(crossinline emBarycenterCoordinates: getHeliocentricEclipticCoordinates,
                                         crossinline moonCoordinates: getHeliocentricEclipticCoordinates): getHeliocentricEclipticCoordinates {
    return { jT ->
        val p = emBarycenterCoordinates(jT)
        val moon = moonCoordinates(jT)
        val earthMoonRatio = 2.7068700387534E7 / 332946.050895
        p - moon * (1.0 / (earthMoonRatio + 1.0))
    }
}

@Heliocentric @Ecliptic @J2000
inline fun createSuspenedEarthMoshierCoordinates(crossinline emBarycenterCoordinates: getHeliocentricEclipticCoordinates,
                                                  crossinline moonCoordinates: getHeliocentricEclipticCoordinates): getCoroutineHeliocentricEclipticCoordinates {
    return { jT ->
        val p = async(CommonPool) { emBarycenterCoordinates(jT) }
        val moon = async(CommonPool) { moonCoordinates(jT) }
        val earthMoonRatio = 2.7068700387534E7 / 332946.050895
         p.await() - moon.await() * (1.0 / (earthMoonRatio + 1.0))
    }
}

fun createLibrationMoshierCoordinates(): (jT: JT) -> Vector {
    return {jT ->
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
        val x = Matrix(Matrix.AXIS_Z, p[2])
        val y = Matrix(Matrix.AXIS_X, p[1])
        val z = Matrix(Matrix.AXIS_Z, p[0])
        val mM = x * y * z

        // Get rotation matrix around x axis with eps
        val mQ = Matrix(Matrix.AXIS_X, 84381.406173 * ARCSEC_TO_RAD)

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
        RectangularVector(p)
    }
}


private fun getMoshierVector(jT: JT, data: MoshierData) =
        gplan(jT, data.args, data.distance, data.tabb, data.tabl,
                data.tabr, data.max_harmonic, data.max_power_of_t,
                data.maxargs, data.timescale, data.trunclvl)