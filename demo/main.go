package main

import (
	"fmt"
	"image/color"
	"log"
	"math"

	"github.com/mowazzem/reedsshepp"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

func plotArrow(p *plot.Plot, x, y, yaw float64, isStart bool) {
	length := 0.5
	arrowHeadSize := 0.2

	x2 := x + length*math.Cos(yaw)
	y2 := y + length*math.Sin(yaw)

	// Draw main arrow line
	line, err := plotter.NewLine(plotter.XYs{{X: x, Y: y}, {X: x2, Y: y2}})
	if err == nil {
		if isStart {
			line.Color = color.RGBA{R: 255, A: 255} // red line
		} else {
			line.Color = color.RGBA{B: 255, A: 255} // blue line
		}
		p.Add(line)
	}

	headAngle := math.Pi / 6 // 30 degrees

	leftAngle := yaw + math.Pi - headAngle
	rightAngle := yaw + math.Pi + headAngle

	xLeft := x2 + arrowHeadSize*math.Cos(leftAngle)
	yLeft := y2 + arrowHeadSize*math.Sin(leftAngle)

	xRight := x2 + arrowHeadSize*math.Cos(rightAngle)
	yRight := y2 + arrowHeadSize*math.Sin(rightAngle)

	leftHead, _ := plotter.NewLine(plotter.XYs{{X: x2, Y: y2}, {X: xLeft, Y: yLeft}})
	rightHead, _ := plotter.NewLine(plotter.XYs{{X: x2, Y: y2}, {X: xRight, Y: yRight}})

	if isStart {
		leftHead.Color = color.RGBA{R: 255, A: 255}
		rightHead.Color = color.RGBA{R: 255, A: 255}
	} else {
		leftHead.Color = color.RGBA{B: 255, A: 255}
		rightHead.Color = color.RGBA{B: 255, A: 255}
	}

	p.Add(leftHead)
	p.Add(rightHead)
}

func main() {
	startX, startY, startYaw := -5.0, -5.0, reedsshepp.Deg2rad(-270)
	endX, endY, endYaw := -5.0, 10.0, reedsshepp.Deg2rad(-180)
	maxCurvature := 2.0
	stepSize := 0.1

	xs, ys, yaws, ctypes, lengths, directions := reedsshepp.ReedsSheppPathPlanning(startX, startY, startYaw, endX, endY, endYaw, maxCurvature, stepSize)
	if xs == nil {
		fmt.Println("Could not find a path.")
		return
	}
	_ = yaws
	_ = ctypes
	_ = lengths

	fmt.Println(endYaw, yaws[len(yaws)-1], reedsshepp.Deg2rad(90)) // this proves the final target and final heading rotation is same

	p := plot.New()
	p.Title.Text = "Reeds-Shepp Path"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"

	// Draw path with forward/backward color
	for i := 0; i < len(xs)-1; i++ {
		lineData := plotter.XYs{
			{X: xs[i], Y: ys[i]},
			{X: xs[i+1], Y: ys[i+1]},
		}
		line, err := plotter.NewLine(lineData)
		if err != nil {
			log.Fatalf("Failed to create path segment: %v", err)
		}
		if directions[i] > 0 {
			line.Color = color.RGBA{G: 200, A: 255} // Green = forward
		} else {
			line.Color = color.RGBA{R: 200, A: 255} // Red = reverse
		}
		line.Width = vg.Points(2)
		p.Add(line)
	}

	plotArrow(p, startX, startY, startYaw, true)
	plotArrow(p, endX, endY, endYaw, false)

	// Dummy lines for legend entries (forward/backward)
	addLegendLine(p, "Forward", color.RGBA{G: 200, A: 255})
	addLegendLine(p, "Backward", color.RGBA{R: 200, A: 255})

	// Save the plot
	if err := p.Save(6*vg.Inch, 6*vg.Inch, "reeds_shepp_path.png"); err != nil {
		log.Fatalf("Failed to save plot: %v", err)
	}

	fmt.Println("Plot saved to 'reeds_shepp_path.png'")
}

func addLegendLine(p *plot.Plot, label string, col color.Color) {
	fn := plotter.NewFunction(func(x float64) float64 { return 0 })
	fn.Color = col
	fn.Width = vg.Points(2)
	p.Legend.Add(label, fn)
}
