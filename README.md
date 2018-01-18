# EnPEn

There it is. The code producing all of my dissertation. It's academic and
therefore it's ugly, won't compile and will lead to people code up the same
stuff, just better (tm).

## What it does

EnPEN solves the coupled Nernst-Planck-Poisson equations in one dimension,
while accounting for chemical reactions. It can handle any ion exchange
membrane setup you can think of (multilayered assemblies, bipolar membranes).
Add some convection, why not.  It can simulate electrical impedance
spectroscopy up to unmeasurable frequencies and down to unmeasurable
frequencies. It's tailored C code, so it'll solve the problem faster and better
than stuff you can buy.

## You think you have a problem that EnPEn can solve?

Don't be afraid to ask. Mail me. I can tell you if EnPEn can do it out of the
box or what parts of the code you might want to change. For instance, change
anything but files in the impl/ directory, or be careful about it.  Since this
is kind of a library, I wrote several python scripts to wrap around it, for
instance: impedance spectroscopy lends itself to having several simulations at
the same time.

## Citing and validation

In the spirit of good scientific practice, you could cite, for steady state
code:

https://www.nature.com/articles/srep11583

and for anything dynamic code:

http://www.sciencedirect.com/science/article/pii/S0376738816311140

## License

Academia? Go ahead. Otherwise ask me.
