package reedsshepp

import (
	"fmt"
	"math"
)

func mod2pi(x float64) float64 {
	v := math.Mod(x, math.Copysign(2.0*math.Pi, x))
	if v < -math.Pi {
		v += 2.0 * math.Pi
	} else {
		v -= -2.0 * math.Pi
	}
	return v
}

func deg2rad(x float64) float64 {
	return x * (math.Pi / 180.0)
}

func Deg2rad(x float64) float64 {
	return deg2rad(x)
}

func rad2deg(x float64) float64 {
	return x * (180.0 * math.Pi)
}

func angleMod(x float64, zero2pi bool, degree bool) float64 {
	if degree {
		x = deg2rad(x)
	}

	var modAngle float64

	if zero2pi {
		modAngle = math.Mod(x, 2*math.Pi)
	} else {
		modAngle = math.Mod((x + math.Pi), (2*math.Pi)-math.Pi)
	}

	if degree {
		modAngle = rad2deg(modAngle)
	}

	return modAngle
}

func pi2Pi(x float64) float64 {
	return angleMod(x, false, false)
}

func sliceAbs(xs []float64) []float64 {
	ys := []float64{}
	for _, x := range xs {
		ys = append(ys, math.Abs(x))
	}
	return ys
}

func sliceSum(xs []float64) float64 {
	ys := 0.0
	for _, x := range xs {
		ys += x
	}
	return ys
}

type Path struct {
	Lengths    []float64 // negative value is backward
	CTypes     []string  // course segment type ,"S": straight, "L": left, "R": right
	L          float64   // total lengths of the path
	X          []float64 // x positions
	Y          []float64 // y positions
	Yaw        []float64 // rotaiton [rad]
	Directions []int     // directions [1: forward, -1: bakckward]
}

func checkTypesEquality(ct1 []string, ct2 []string) bool {
	mp := map[string]struct{}{}
	for _, ct := range ct1 {
		mp[ct] = struct{}{}
	}
	for _, ct := range ct2 {
		if _, ok := mp[ct]; !ok {
			return false
		}
	}
	return true
}

func checkLengthsAreLess(lengths []float64, total float64, stepSize float64) bool {
	y := sliceSum(sliceAbs(lengths))
	return (y - total) <= stepSize
}

func setPath(paths []Path, lengths []float64, ctypes []string, stepSize float64) []Path {
	path := Path{
		CTypes:  ctypes,
		Lengths: lengths,
		L:       sliceSum(sliceAbs(lengths)),
	}
	for _, iPath := range paths {
		isTypeSame := checkTypesEquality(iPath.CTypes, path.CTypes)
		isLengthClose := checkLengthsAreLess(iPath.Lengths, path.L, stepSize)
		if isTypeSame && isLengthClose {
			return paths
		}
	}
	paths = append(paths, path)
	return paths
}

func polar(x, y float64) (float64, float64) {
	r := math.Hypot(x, y)
	theta := math.Atan2(y, x)
	return r, theta
}

func leftStraightLeft(x, y, phi float64) (bool, []float64, []string) {
	dx := x - math.Sin(phi)
	dy := y - 1.0 + math.Cos(phi)

	u, t := polar(dx, dy)

	if t >= 0.0 && t <= math.Pi {
		v := mod2pi(phi - t)
		if v >= 0.0 && v <= math.Pi {
			return true, []float64{t, u, v}, []string{"L", "S", "L"}
		}
	}
	return false, nil, nil
}

func leftStraightRight(x, y, phi float64) (bool, []float64, []string) {
	dx := x + math.Sin(phi)
	dy := y - 1.0 - math.Cos(phi)

	u1, t1 := polar(dx, dy)
	u1 = u1 * u1 // square of u1

	if u1 >= 4.0 {
		u := math.Sqrt(u1 - 4.0)
		theta := math.Atan2(2.0, u)
		t := mod2pi(t1 + theta)
		v := mod2pi(t - phi)

		if t >= 0.0 && v >= 0.0 {
			return true, []float64{t, u, v}, []string{"L", "S", "R"}
		}
	}

	return false, nil, nil
}

func leftXRightXLeft(x, y, phi float64) (bool, []float64, []string) {
	zeta := x - math.Sin(phi)
	eeta := y - 1.0 + math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 <= 4.0 {
		A := math.Acos(0.25 * u1)
		t := mod2pi(A + theta + math.Pi/2)
		u := mod2pi(math.Pi - 2.0*A)
		v := mod2pi(phi - t - u)
		return true, []float64{t, -u, v}, []string{"L", "R", "L"}
	}

	return false, nil, nil
}

func leftXRightLeft(x, y, phi float64) (bool, []float64, []string) {
	zeta := x - math.Sin(phi)
	eeta := y - 1.0 + math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 <= 4.0 {
		A := math.Acos(0.25 * u1)
		t := mod2pi(A + theta + math.Pi/2)
		u := mod2pi(math.Pi - 2.0*A)
		v := mod2pi(-phi + t + u)
		return true, []float64{t, -u, -v}, []string{"L", "R", "L"}
	}

	return false, nil, nil
}

func leftRightXLeft(x, y, phi float64) (bool, []float64, []string) {
	zeta := x - math.Sin(phi)
	eeta := y - 1.0 + math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 <= 4.0 {
		u := math.Acos(1.0 - 0.125*u1*u1)
		A := math.Asin(2.0 * math.Sin(u) / u1)
		t := mod2pi(-A + theta + math.Pi/2)
		v := mod2pi(t - u - phi)
		return true, []float64{t, u, -v}, []string{"L", "R", "L"}
	}

	return false, nil, nil
}

func leftRightXLeftRight(x, y, phi float64) (bool, []float64, []string) {
	zeta := x + math.Sin(phi)
	eeta := y - 1.0 - math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	// Solutions with 2 < u1 <= 4 are considered sub-optimal
	// Solutions do not exist for u1 > 4
	if u1 <= 2.0 {
		A := math.Acos((u1 + 2.0) * 0.25)
		t := mod2pi(theta + A + math.Pi/2)
		u := mod2pi(A)
		v := mod2pi(phi - t + 2.0*u)

		if t >= 0.0 && u >= 0.0 && v >= 0.0 {
			return true, []float64{t, u, -u, -v}, []string{"L", "R", "L", "R"}
		}
	}

	return false, nil, nil
}

func leftXRightLeftXRight(x, y, phi float64) (bool, []float64, []string) {
	zeta := x + math.Sin(phi)
	eeta := y - 1.0 - math.Cos(phi)

	u1, theta := polar(zeta, eeta)
	u2 := (20.0 - u1*u1) / 16.0

	if u2 >= 0.0 && u2 <= 1.0 {
		u := math.Acos(u2)
		A := math.Asin(2.0 * math.Sin(u) / u1)
		t := mod2pi(theta + A + math.Pi/2)
		v := mod2pi(t - phi)

		if t >= 0.0 && v >= 0.0 {
			return true, []float64{t, -u, -u, v}, []string{"L", "R", "L", "R"}
		}
	}

	return false, nil, nil
}

func leftXRight90StraightLeft(x, y, phi float64) (bool, []float64, []string) {
	zeta := x - math.Sin(phi)
	eeta := y - 1.0 + math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 >= 2.0 {
		u := math.Sqrt(u1*u1-4) - 2
		A := math.Atan2(2, math.Sqrt(u1*u1-4))
		t := mod2pi(theta + A + math.Pi/2)
		v := mod2pi(t - phi + math.Pi/2)
		if t >= 0 && v >= 0 {
			return true, []float64{t, -math.Pi / 2, -u, -v}, []string{"L", "R", "S", "L"}
		}
	}

	return false, nil, nil
}

func leftStraightRight90XLeft(x, y, phi float64) (bool, []float64, []string) {
	zeta := x - math.Sin(phi)
	eeta := y - 1.0 + math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 >= 2.0 {
		u := math.Sqrt(u1*u1-4.0) - 2.0
		A := math.Atan2(math.Sqrt(u1*u1-4.0), 2.0)
		t := mod2pi(theta - A + math.Pi/2)
		v := mod2pi(t - phi - math.Pi/2)
		if t >= 0 && v >= 0 {
			return true, []float64{t, u, math.Pi / 2, -v}, []string{"L", "S", "R", "L"}
		}
	}

	return false, nil, nil
}

func leftXRight90StraightRight(x, y, phi float64) (bool, []float64, []string) {
	zeta := x + math.Sin(phi)
	eeta := y - 1.0 - math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 >= 2.0 {
		t := mod2pi(theta + math.Pi/2)
		u := u1 - 2.0
		v := mod2pi(phi - t - math.Pi/2)
		if t >= 0 && v >= 0 {
			return true, []float64{t, -math.Pi / 2, -u, -v}, []string{"L", "R", "S", "R"}
		}
	}

	return false, nil, nil
}

func leftStraightLeft90XRight(x, y, phi float64) (bool, []float64, []string) {
	zeta := x + math.Sin(phi)
	eeta := y - 1.0 - math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 >= 2.0 {
		t := mod2pi(theta)
		u := u1 - 2.0
		v := mod2pi(phi - t - math.Pi/2)
		if t >= 0 && v >= 0 {
			return true, []float64{t, u, math.Pi / 2, -v}, []string{"L", "S", "L", "R"}
		}
	}

	return false, nil, nil
}

func leftXRight90StraightLeft90XRight(x, y, phi float64) (bool, []float64, []string) {
	zeta := x + math.Sin(phi)
	eeta := y - 1.0 - math.Cos(phi)

	u1, theta := polar(zeta, eeta)

	if u1 >= 4.0 {
		u := math.Sqrt(u1*u1-4) - 4
		A := math.Atan2(2, math.Sqrt(u1*u1-4))
		t := mod2pi(theta + A + math.Pi/2)
		v := mod2pi(t - phi)

		if t >= 0 && v >= 0 {
			return true, []float64{t, -math.Pi / 2, -u, -math.Pi / 2, v}, []string{"L", "R", "S", "L", "R"}
		}
	}

	return false, nil, nil
}

func timeflip(travelDistances []float64) []float64 {
	flipped := make([]float64, len(travelDistances))
	for i, v := range travelDistances {
		flipped[i] = -v
	}
	return flipped
}

func reflect(steeringDirections []string) []string {
	switchDir := func(dir string) string {
		switch dir {
		case "L":
			return "R"
		case "R":
			return "L"
		default:
			return "S"
		}
	}

	reflected := make([]string, len(steeringDirections))
	for i, dir := range steeringDirections {
		reflected[i] = switchDir(dir)
	}
	return reflected
}

func generatePath(q0, q1 [3]float64, maxCurvature, stepSize float64) []Path {
	dx := q1[0] - q0[0]
	dy := q1[1] - q0[1]
	dth := q1[2] - q0[2]
	c := math.Cos(q0[2])
	s := math.Sin(q0[2])
	x := (c*dx + s*dy) * maxCurvature
	y := (-s*dx + c*dy) * maxCurvature
	stepSize *= maxCurvature

	var paths []Path

	// List of path functions with signature:
	// func(float64, float64, float64) (bool, []float64, []string)
	pathFunctions := []func(float64, float64, float64) (bool, []float64, []string){
		leftStraightLeft, leftStraightRight,
		leftXRightXLeft, leftXRightLeft, leftRightXLeft,
		leftRightXLeftRight, leftXRightLeftXRight,
		leftXRight90StraightLeft, leftXRight90StraightRight,
		leftStraightRight90XLeft, leftStraightLeft90XRight,
		leftXRight90StraightLeft90XRight,
	}

	for _, pathFunc := range pathFunctions {
		flag, travelDistances, steeringDirs := pathFunc(x, y, dth)
		if flag {
			sumAbs := 0.0
			for _, d := range travelDistances {
				sumAbs += math.Abs(d)
			}
			for _, distance := range travelDistances {
				if 0.1*sumAbs < math.Abs(distance) && math.Abs(distance) < stepSize {
					fmt.Println("Step size too large for Reeds-Shepp paths.")
					return nil
				}
			}
			paths = setPath(paths, travelDistances, steeringDirs, stepSize)
		}

		flag, travelDistances, steeringDirs = pathFunc(-x, y, -dth)
		if flag {
			sumAbs := 0.0
			for _, d := range travelDistances {
				sumAbs += math.Abs(d)
			}
			for _, distance := range travelDistances {
				if 0.1*sumAbs < math.Abs(distance) && math.Abs(distance) < stepSize {
					fmt.Println("Step size too large for Reeds-Shepp paths.")
					return nil
				}
			}
			travelDistances = timeflip(travelDistances)
			paths = setPath(paths, travelDistances, steeringDirs, stepSize)
		}

		flag, travelDistances, steeringDirs = pathFunc(x, -y, -dth)
		if flag {
			sumAbs := 0.0
			for _, d := range travelDistances {
				sumAbs += math.Abs(d)
			}
			for _, distance := range travelDistances {
				if 0.1*sumAbs < math.Abs(distance) && math.Abs(distance) < stepSize {
					fmt.Println("Step size too large for Reeds-Shepp paths.")
					return nil
				}
			}
			steeringDirs = reflect(steeringDirs)
			paths = setPath(paths, travelDistances, steeringDirs, stepSize)
		}

		flag, travelDistances, steeringDirs = pathFunc(-x, -y, dth)
		if flag {
			sumAbs := 0.0
			for _, d := range travelDistances {
				sumAbs += math.Abs(d)
			}
			for _, distance := range travelDistances {
				if 0.1*sumAbs < math.Abs(distance) && math.Abs(distance) < stepSize {
					fmt.Println("Step size too large for Reeds-Shepp paths.")
					return nil
				}
			}
			travelDistances = timeflip(travelDistances)
			steeringDirs = reflect(steeringDirs)
			paths = setPath(paths, travelDistances, steeringDirs, stepSize)
		}
	}

	return paths
}

func calcInterpolateDistsList(lengths []float64, stepSize float64) [][]float64 {
	interpolateDistsList := make([][]float64, 0, len(lengths))

	for _, length := range lengths {
		var dDist float64
		if length >= 0 {
			dDist = stepSize
		} else {
			dDist = -stepSize
		}

		var interpDists []float64
		for dist := 0.0; (dDist > 0 && dist < length) || (dDist < 0 && dist > length); dist += dDist {
			interpDists = append(interpDists, dist)
		}
		interpDists = append(interpDists, length)

		interpolateDistsList = append(interpolateDistsList, interpDists)
	}

	return interpolateDistsList
}

func interpolate(dist, length float64, mode string, maxCurvature float64, originX, originY, originYaw float64) (float64, float64, float64, int) {
	var x, y, yaw float64

	if mode == "S" {
		x = originX + dist/maxCurvature*math.Cos(originYaw)
		y = originY + dist/maxCurvature*math.Sin(originYaw)
		yaw = originYaw
	} else {
		ldx := math.Sin(dist) / maxCurvature
		ldy := 0.0

		if mode == "L" {
			ldy = (1.0 - math.Cos(dist)) / maxCurvature
			yaw = originYaw + dist
		} else if mode == "R" {
			ldy = (1.0 - math.Cos(dist)) / -maxCurvature
			yaw = originYaw - dist
		}

		gdx := math.Cos(-originYaw)*ldx + math.Sin(-originYaw)*ldy
		gdy := -math.Sin(-originYaw)*ldx + math.Cos(-originYaw)*ldy

		x = originX + gdx
		y = originY + gdy
	}

	var direction int
	if length > 0 {
		direction = 1
	} else {
		direction = -1
	}

	return x, y, yaw, direction
}

func generateLocalCourse(lengths []float64, modes []string, maxCurvature, stepSize float64) ([]float64, []float64, []float64, []int) {
	interpolateDistsList := calcInterpolateDistsList(lengths, stepSize*maxCurvature)

	originX, originY, originYaw := 0.0, 0.0, 0.0

	var xs, ys, yaws []float64
	var directions []int

	for i, interpDists := range interpolateDistsList {
		mode := modes[i]
		length := lengths[i]

		for _, dist := range interpDists {
			x, y, yaw, direction := interpolate(dist, length, mode, maxCurvature, originX, originY, originYaw)
			xs = append(xs, x)
			ys = append(ys, y)
			yaws = append(yaws, yaw)
			directions = append(directions, direction)
		}

		originX = xs[len(xs)-1]
		originY = ys[len(ys)-1]
		originYaw = yaws[len(yaws)-1]
	}

	return xs, ys, yaws, directions
}

func calcPaths(sx, sy, syaw, gx, gy, gyaw, maxCurvature, stepSize float64) []Path {
	q0 := [3]float64{sx, sy, syaw}
	q1 := [3]float64{gx, gy, gyaw}

	paths := generatePath(q0, q1, maxCurvature, stepSize)

	for i := range paths {
		path := &paths[i]

		xs, ys, yaws, directions := generateLocalCourse(path.Lengths, path.CTypes, maxCurvature, stepSize)

		// convert global coordinate
		cosTheta := math.Cos(-q0[2])
		sinTheta := math.Sin(-q0[2])

		path.X = make([]float64, len(xs))
		path.Y = make([]float64, len(ys))
		path.Yaw = make([]float64, len(yaws))

		for j := range xs {
			ix, iy := xs[j], ys[j]
			path.X[j] = cosTheta*ix + sinTheta*iy + q0[0]
			path.Y[j] = -sinTheta*ix + cosTheta*iy + q0[1]
			path.Yaw[j] = pi2Pi(yaws[j] + q0[2])
		}

		path.Directions = directions

		for j := range path.Lengths {
			path.Lengths[j] /= maxCurvature
		}
		path.L /= maxCurvature
	}

	return paths
}

func ReedsSheppPathPlanning(sx, sy, syaw, gx, gy, gyaw, maxCurvature, stepSize float64) ([]float64, []float64, []float64, []string, []float64, []int) {
	paths := calcPaths(sx, sy, syaw, gx, gy, gyaw, maxCurvature, stepSize)
	if len(paths) == 0 {
		return nil, nil, nil, nil, nil, nil // could not generate any path
	}

	bestIndex := 0
	minCost := math.Abs(paths[0].L)
	for i, p := range paths {
		if cost := math.Abs(p.L); cost < minCost {
			minCost = cost
			bestIndex = i
		}
	}

	bPath := paths[bestIndex]

	return bPath.X, bPath.Y, bPath.Yaw, bPath.CTypes, bPath.Lengths, bPath.Directions
}
