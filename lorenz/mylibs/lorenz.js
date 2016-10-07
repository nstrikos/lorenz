function f0(x, y0, y1, y2, sigma)
{
    return sigma * (y1 - y0);
}

function f1(x, y0, y1, y2, rho)
{
    return y0 * (rho - y2) - y1;
}

function f2(x, y0, y1, y2, beta)
{
    return y0 * y1 - beta * y2;
}



function rk4(x, xstop, h, y, a, aRows)
{

    /*
    x -             start  time
        xstop -         final time
        h -             time step
        y -             initial values
        y0, y1, y2 -    m1 initial position
        y3, y4, y5 -    m2 initial position
        y6, y7, y8 -    m1 initial velocity
        y9, y10, y11 -  m1 initial velocity
        m1 -            m1 mass
        m2 -            m2 mass
        G -             gravitational constant
        a -             array of solutions
    */

    /* k matrix */

    var sigma = 5.0;
    var rho = 15.0;
    var beta = 1.0;

    var i = 0;
    var k1_1 = 0;
    var k1_2 = 0;
    var k1_3 = 0;

    var k2_1 = 0;
    var k2_2 = 0;
    var k2_3 = 0;

    var k3_1 = 0;
    var k3_2 = 0;
    var k3_3 = 0;

    var k4_1 = 0;
    var k4_2 = 0;
    var k4_3 = 0;

    for (i=0; i<aRows; i++)
    {
        if (xstop - x < h)
            h = xstop - x;

        k1_1 = h * f0(x, y[0], y[1], y[2], sigma);
        k1_2 = h * f1(x, y[0], y[1], y[2], rho);
        k1_3 = h * f2(x, y[0], y[1], y[2], beta);

        k2_1 = h * f0(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, sigma);
        k2_2 = h * f1(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, rho);
        k2_3 = h * f2(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, beta);

        k3_1 = h * f0(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, sigma);
        k3_2 = h * f1(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, rho);
        k3_3 = h * f2(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, beta);

        k4_1 = h * f0(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, sigma);
        k4_2 = h * f1(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, rho);
        k4_3 = h * f2(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, beta);

        y[0] += (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
        y[1] += (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
        y[2] += (k1_3 + 2*k2_3 + 2*k3_3 + k4_3)/6;

        a[i][0] = y[0];
        a[i][1] = y[1];
        a[i][2] = y[2];

        x += h;
    }
}
	
function rk4_b(y, h, a )
	{
	
		/*
			y -             array of initial values
			h -             time step
			a -             array of solutions
		*/

		/* k matrix */

		var sigma = 5.0;
		var rho = 15.0;
		var beta = 1.0;

		var i = 0;
		var k1_1 = 0;
		var k1_2 = 0;
		var k1_3 = 0;

		var k2_1 = 0;
		var k2_2 = 0;
		var k2_3 = 0;

		var k3_1 = 0;
		var k3_2 = 0;
		var k3_3 = 0;

		var k4_1 = 0;
		var k4_2 = 0;
		var k4_3 = 0;

		
		k1_1 = h * f0(x, y[0], y[1], y[2], sigma);
		k1_2 = h * f1(x, y[0], y[1], y[2], rho);
		k1_3 = h * f2(x, y[0], y[1], y[2], beta);

		k2_1 = h * f0(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, sigma);
		k2_2 = h * f1(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, rho);
		k2_3 = h * f2(x + h/2, y[0] + k1_1/2, y[1] + k1_2/2, y[2] + k1_3/2, beta);

		k3_1 = h * f0(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, sigma);
		k3_2 = h * f1(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, rho);
		k3_3 = h * f2(x + h/2, y[0] + k2_1/2, y[1] + k2_2/2, y[2] + k2_3/2, beta);

		k4_1 = h * f0(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, sigma);
		k4_2 = h * f1(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, rho);
		k4_3 = h * f2(x + h, y[0] + k3_1, y[1] + k3_2, y[2] + k3_3, beta);

		y[0] += (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
		y[1] += (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
		y[2] += (k1_3 + 2*k2_3 + 2*k3_3 + k4_3)/6;

		a[0] = y[0];
		a[1] = y[1];
		a[2] = y[2];					
	}

