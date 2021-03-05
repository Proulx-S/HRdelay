x = randn(30);
imagesc(x)
colorbar
[COEFF1, SCORE1, LATENT1, TSQUARED1, EXPLAINED1, MU1] = pca(x);
