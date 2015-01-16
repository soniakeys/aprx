// Public domain.

// Implementation of "Keplerian Elements for Approximate Positions of the
// Major Planets" by  E M Standish.
//
// Approximate errors, from http://ssd.jpl.nasa.gov/?planet_pos
//
//                     1800 -- 2050             3000 BCE to 3000 CE
//              ---------------------------  --------------------------
//                 RA      Dec.       r         RA      Dec.       r
//              (arcsec) (arcsec) (1000 km)	 (arcsec) (arcsec) (1000 km)
//  Mercury        15        1        1         20       15        1
//  Venus          20        1        4         40       30        8
//  Earth-Moon
//  Barycenter     20        8        6         40       15       15
//  Mars           40        2       25        100       40       30
//  Jupiter       400       10      600        600      100     1000
//  Saturn        600       25     1500       1000      100     4000
//  Uranus         50        2     1000       2000       30     8000
//  Neptune        10        1      200        400       15     4000
//  Pluto           5        2      300        400      100     2500
package aprx

// aprx_pos_planets.pdf, p_elem_t1.txt, p_elem_t2.txt accessed 8 Feb 2012
// from http://ssd.jpl.nasa.gov/?planet_pos

import (
	"errors"
	"math"
	"strconv"
	"strings"
)

// Planet constants
const (
	Mercury = iota
	Venus
	EMBary
	Mars
	Jupiter
	Saturn
	Uranus
	Neptune
	Pluto
	nPlanets // do not change
)

// Julian dates corresponding to years
const (
	j3000b = 625673.5
	j1800  = 2378496.5
	j2000  = 2451545.0
	j2050  = 2469807.5
	j3000  = 2816787.5
)

// EclPos returns approximate ecliptic coordinates.
//
// Parameter planet is one of the planet constants.
// Parameter tEph is a Julian date in the range 2378496.5 to 2469807.5,
// which corresponds to the years 1800 to 2050 CE.
//
// Returned cartesian coordinates are in AU.
func EclPos(planet int, tEph float64) (x, y, z float64, err error) {
	if planet < 0 || planet >= nPlanets {
		err = errors.New("Invalid planet.")
		return
	}
	if tEph < j1800 || tEph > j2050 {
		err = errors.New("Date out of valid interval.")
		return
	}
	tCen := (tEph - j2000) / 36525
	k := rates[planet].at(tCen)
	x, y, z = eclM(k, k.l-k.ϖ)
	return
}

// EclPosLong returns approximate ecliptic coordinates using "long interval."
//
// Parameter planet is one of the planet constants.
// Parameter tEph is a Julian date in the range 3000 BCE to 3000 CE.
//
// Returned cartesian coordinates are in AU.
func EclPosLong(planet int, tEph float64) (x, y, z float64, err error) {
	if planet < 0 || planet >= nPlanets {
		err = errors.New("Invalid planet.")
		return
	}
	if tEph < j3000b || tEph > j3000 {
		err = errors.New("Date out of valid interval.")
		return
	}
	tCen := (tEph - j2000) / 36525
	k := longRates[planet].at(tCen)
	ma := k.l - k.ϖ
	if planet >= Jupiter {
		px := planet - Jupiter
		ma += b[px] * tCen * tCen
		if planet < Pluto {
			sin, cos := math.Sincos(f[px] * tCen * math.Pi / 180)
			ma += s[px]*sin + c[px]*cos
		}
	}
	x, y, z = eclM(k, ma)
	return
}

// obliquity of ecliptic at J2000
var sε, cε = math.Sincos(23.43928 * math.Pi / 180)

// Equ returns equatorial coordinates in the J2000 frame, given
// ecliptic coordinates.  Note that it has no way to check that
// arguments really represent ecliptic coordinates.
func Equ(x, y, z float64) (float64, float64, float64) {
	return x, cε*y - sε*z, sε*y + cε*z
}

func eclM(k *kep, ma float64) (x, y, z float64) {
	se, ce := math.Sincos(k.ea(ma) * math.Pi / 180)
	xp := k.a * (ce - k.e)
	yp := k.a * math.Sqrt(1-k.e*k.e) * se
	sω, cω := math.Sincos(k.ω() * math.Pi / 180)
	sΩ, cΩ := math.Sincos(k.Ω * math.Pi / 180)
	si, ci := math.Sincos(k.i * math.Pi / 180)
	x = (cω*cΩ-sω*sΩ*ci)*xp + (-sω*cΩ-cω*sΩ*ci)*yp
	y = (cω*sΩ+sω*cΩ*ci)*xp + (-sω*sΩ+cω*cΩ*ci)*yp
	z = sω*si*xp + cω*si*yp
	return
}

type eleRate struct {
	ele, rate float64
}

// compute element at time t, where t is centuries from J2000.
func (e *eleRate) at(tCen float64) float64 {
	return e.ele + e.rate*tCen
}

type kepRates struct {
	a eleRate // semi-major axis in au, au/day
	e eleRate // eccentricity, unitless, /day
	i eleRate // inclination in degrees, degrees/day
	l eleRate // mean longitude in degrees, degrees/day
	ϖ eleRate // longitude of perihelion (ω+Ω) in degrees, degrees/day
	Ω eleRate // longitude of ascending node in degrees, degrees/day
}

func (k *kepRates) at(t float64) *kep {
	return &kep{
		a: k.a.at(t),
		e: k.e.at(t),
		i: k.i.at(t),
		l: k.l.at(t),
		ϖ: k.ϖ.at(t),
		Ω: k.Ω.at(t),
	}
}

type kep struct {
	a float64 // semi-major axis in au
	e float64 // eccentricity, unitless
	i float64 // inclination in degrees
	l float64 // mean longitude in degrees
	ϖ float64 // longitude of perihelion (ω+Ω) in degrees
	Ω float64 // longitude of ascending node in degrees
}

// argument of perihelion
func (k *kep) ω() float64 {
	return k.ϖ - k.Ω
}

// eccentric anomaly
func (k *kep) ea(ma float64) float64 {
	m := math.Mod(ma, 360)
	if m < -180 {
		m += 360
	} else if m > 180 {
		m -= 360
	}
	const tol = 1e-6
	e0 := m - k.e*math.Sin(m*math.Pi/180)
	eRad := k.e * math.Pi / 180
	for {
		se, ce := math.Sincos(e0 * math.Pi / 180)
		dm := m - (e0 - k.e*se)
		de := dm / (1 - eRad*ce)
		e0 += de
		if math.Abs(de) < tol {
			break
		}
	}
	return e0
}

var rates, longRates [nPlanets]kepRates

func init() {
	// data copied from Table 1, file p_elem_t1.txt.  A table of Go numeric
	// literals would be convenient, but data is copied as a single large
	// string to avoid typographical errors.
	//         a              e               I                L            long.peri.      long.node.
	//     AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
	parseRates(&rates, `
Mercury   0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593
          0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081
Venus     0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255
          0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418
EM Bary   1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0
          0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0
Mars      1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
          0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343
Jupiter   5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909
         -0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106
Saturn    9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448
         -0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794
Uranus   19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503
         -0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589
Neptune  30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574
          0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664
Pluto    39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684
         -0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482
`)
	// data copied from Table 2a, p_elem_t2.txt
	//         a              e               I                L            long.peri.      long.node.
	//     AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy     deg, deg/Cy
	parseRates(&longRates, `
Mercury   0.38709843      0.20563661      7.00559432      252.25166724     77.45771895     48.33961819
          0.00000000      0.00002123     -0.00590158   149472.67486623      0.15940013     -0.12214182
Venus     0.72332102      0.00676399      3.39777545      181.97970850    131.76755713     76.67261496
         -0.00000026     -0.00005107      0.00043494    58517.81560260      0.05679648     -0.27274174
EM Bary   1.00000018      0.01673163     -0.00054346      100.46691572    102.93005885     -5.11260389
         -0.00000003     -0.00003661     -0.01337178    35999.37306329      0.31795260     -0.24123856
Mars      1.52371243      0.09336511      1.85181869       -4.56813164    -23.91744784     49.71320984
          0.00000097      0.00009149     -0.00724757    19140.29934243      0.45223625     -0.26852431
Jupiter   5.20248019      0.04853590      1.29861416       34.33479152     14.27495244    100.29282654
         -0.00002864      0.00018026     -0.00322699     3034.90371757      0.18199196      0.13024619
Saturn    9.54149883      0.05550825      2.49424102       50.07571329     92.86136063    113.63998702
         -0.00003065     -0.00032044      0.00451969     1222.11494724      0.54179478     -0.25015002
Uranus   19.18797948      0.04685740      0.77298127      314.20276625    172.43404441     73.96250215
         -0.00020455     -0.00001550     -0.00180155      428.49512595      0.09266985      0.05739699
Neptune  30.06952752      0.00895439      1.77005520      304.22289287     46.68158724    131.78635853
          0.00006447      0.00000818      0.00022400      218.46515314      0.01009938     -0.00606302
Pluto    39.48686035      0.24885238     17.14104260      238.96535011    224.09702598    110.30167986
          0.00449751      0.00006016      0.00000501      145.18042903     -0.00968827     -0.00809981
`)
	// data copied from Table 2b, p_elem_t2.txt
	//          b             c             s            f
	lines := strings.Split(`
Jupiter   -0.00012452    0.06064060   -0.35635438   38.35125000
Saturn     0.00025899   -0.13434469    0.87320147   38.35125000
Uranus     0.00058331   -0.97731848    0.17689245    7.67025000
Neptune   -0.00041348    0.68346318   -0.10162547    7.67025000
Pluto     -0.01262724`, "\n")[1:]

	for p, line := range lines {
		fld := strings.Fields(line)
		x, err := strconv.ParseFloat(fld[1], 64)
		if err != nil {
			panic(err)
		}
		b[p] = x
		if p+Jupiter == Pluto {
			return
		}
		if x, err = strconv.ParseFloat(fld[2], 64); err != nil {
			panic(err)
		}
		c[p] = x
		if x, err = strconv.ParseFloat(fld[3], 64); err != nil {
			panic(err)
		}
		s[p] = x
		if x, err = strconv.ParseFloat(fld[4], 64); err != nil {
			panic(err)
		}
		f[p] = x
	}
}

var (
	b [nPlanets - Jupiter]float64
	c [Pluto - Jupiter]float64
	s [Pluto - Jupiter]float64
	f [Pluto - Jupiter]float64
)

func parseRates(r *[nPlanets]kepRates, table string) {
	lines := strings.Split(table, "\n")[1:]
	for pl := range r {
		sEle := strings.Fields(lines[pl*2][8:]) // 8: skips planet name
		sRate := strings.Fields(lines[pl*2+1][8:])
		r[pl].a = parseRate(sEle[0], sRate[0])
		r[pl].e = parseRate(sEle[1], sRate[1])
		r[pl].i = parseRate(sEle[2], sRate[2])
		r[pl].l = parseRate(sEle[3], sRate[3])
		r[pl].ϖ = parseRate(sEle[4], sRate[4])
		r[pl].Ω = parseRate(sEle[5], sRate[5])
	}
}

func parseRate(sEle, sRate string) (r eleRate) {
	var err error
	if r.ele, err = strconv.ParseFloat(sEle, 64); err != nil {
		panic(err)
	}
	if r.rate, err = strconv.ParseFloat(sRate, 64); err != nil {
		panic(err)
	}
	return
}
