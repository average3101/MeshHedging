
# MeshHedging

**MeshHedging** is a Java library for computing (approximately) optimal policies for portfolio hedging problems, using stochastic mesh methods.

It was developped in the context of my PhD studies, under the supervision of [Pierre L'Ecuyer](http://www-labs.iro.umontreal.ca/~lecuyer/) 
at the [Department of Computer Science and Operations Research](http://en.diro.umontreal.ca/home/) of Université de Montréal. 

## Documentation 

A short presentation of the classes of the MeshHedging library is available in section 5.1
of my doctoral [thesis](https://github.com/average3101/MeshHedging/blob/master/PATremblay_thesis_meshhedging.pdf).


## Author

* **Pierre-Alexandre Tremblay**


## License

This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details

## Dependencies
**MeshHedging** depends on the following external librairies:

* [SSJ](http://simul.iro.umontreal.ca/ssj/)
* [COLT](http://dst.lbl.gov/ACSSoftware/colt/)
* [Commons Math](http://commons.apache.org/proper/commons-math/)

The [Fmin](https://www1.fpl.fs.fed.us/Fmin.java) class by Steve Verrill is also used, but has been integrated in the code of *BrentSolver* class of **MeshHedging**.



