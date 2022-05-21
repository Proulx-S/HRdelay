function theta = wrapToPi(theta)
theta = mod(theta+pi, 2*pi)-pi;