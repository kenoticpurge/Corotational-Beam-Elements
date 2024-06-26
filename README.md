# Corotational-Beam-Elements
You ever just wanted to do some large deformation, nonlinear analysis of beams?

Well you're in luck

This is just some code for beam elements which use the "Corotational Formulation"

You'll find both planar and spatial beam elements, statics and dynamics, and a few cases studies like Canitlever Beams or Rotating Beams

The "Corotational Formulation" used can be found in the following key papers:
  "Efﬁcient formulation for dynamics of corotational 2D beams", "Co-rotational beam elements with warping eﬀects in instability problems", "Consistent co-rotational framework for Euler-Bernoulli and Timoshenko beam-column elements under distributed member loads"
  
The implementation is verified by results in ABAQUS and comparing with the papers.

Use the code however you like, hopefully someone will find it useful.

My advice is to go through the Cantilever Beam folder and run the script "Cantilever_Beam_HHT". Do this while reading through the equations in the paper "Efﬁcient formulation for dynamics of corotational 2D beams"

Large Deformation Cantilever Beam Results:
![untitled](https://github.com/kenoticpurge/Corotational-Beam-Elements/assets/157754800/85a9475a-375b-4190-81ff-f868a5edab5b)

