//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "QMCMikesBracketingStepLengthSelector.h"

double QMCMikesBracketingStepLengthSelector::stepLength(
    QMCObjectiveFunction *function, Array1D<double> & position,
    Array1D<double> & searchDirection, Array1D<double> & gradient, 
    double functionValue)

{
  double initial_step_length=0.001;
  double initial_alpha=0.0;

  //bracket a min in the direction of p starting at position
  bracket(function, position,searchDirection,initial_alpha,
	  initial_step_length,0);

  print_interval();

  //this is the sign that we have recursed too many times 
  //and alpha has gotten so small we are giving up
  //alpha=0.0 is the flag for this to stop
  if(alpha_middle<0.0)
    {
      return 0.0;
    }

  //while we are not done do the following
  bool done=false;
  int counter=0;
  while(!done)
    {
      //rebracket on the previous bracketed min to make it tighter
      quadratic_rebracketer(function,position,searchDirection);

      print_interval();

      //decide if we should terminate this loop
      if(counter>=3)
        {
          done=true;
        }
      counter++;
    }
  return alpha_middle; //the best alpha thus far
}

void QMCMikesBracketingStepLengthSelector::print_interval()
{
  cout<<endl;
  cout<<"alphas:\t"<<alpha_left<<"\t"<<alpha_middle<<"\t"<<alpha_right<<endl;
  cout<<"scores:\t"<<score_left<<"\t"<<score_middle<<"\t"<<score_right<<endl;
  cout<<endl;
}

//This will make the alphas properly bracket a minimum.
void QMCMikesBracketingStepLengthSelector::bracket(QMCObjectiveFunction *OF, 
    Array1D<double> & position,Array1D<double> & searchDirection,
    double alpha_zero,double scale,int recursion_depth)
{
  cout<<"\nbracket()"<<endl;
  cout<<"x: "<<position;
  cout<<"p: "<<searchDirection;
  cout<<alpha_zero<<"\t"<<scale<<"\t"<<recursion_depth<<endl;

  if(recursion_depth<4)
    {
      //how many steps will we try to bracket
      int steps=8;
      
      //allocate a population to be evaulated
      Array1D < Array1D <double> > population;
      population.allocate(steps);
      for(int i=0;i<steps;i++)
        {
          population(i).allocate(position.dim1());
        }
      
      //allocate the alphas to be tried
      Array1D<double> alphas;
      alphas.allocate(steps);
      
      //allocate the scores to be tried
      Array1D<double> scores;
      scores.allocate(steps);
      
      //allocate the results from an objective function evaulation
      Array1D< QMCObjectiveFunctionResult > results;
      results.allocate(steps);

      //determine what alphas we will try
      alphas(0)=alpha_zero;
      double temp=1.0;  //will keep track of 2^n factor
      for(int i=1;i<population.dim1();i++)
        {  
          //do alphas that are 2^n times the search direction
          temp=temp*2.0;
          alphas(i)=alpha_zero+scale*temp;
        }
      
      //fill the population with these trial steps
      for(int i=0;i<population.dim1();i++)
        {
          for(int j=0;j<population(i).dim1();j++)
            {
              population(i)(j)=position(j)+alphas(i)*
                searchDirection(j);
            }
          cout<<population(i);
        }

      cout<<"evaluating in bracketer"<<endl;
      results = OF->evaluate(population);
      cout<<"just after evaluating in bracketer"<<endl;

      //put the scores into the scores vector
      for(int i=0;i<population.dim1();i++)
        {
          scores(i)=results(i).getScore();

          if(IeeeMath::isnan(scores(i)) != 0) 
            {
              scores(i) = 1e30;
            }
        }

      cout << "\nALPHAS = " << alphas;
      cout << "SCORES = " << scores<<endl;;
      
      //find the minimum in the population
      int min_index=0;
      double min_score=scores(0)+1e6;
      for(int i=0;i<population.dim1();i++)
        {
          if(min_score>scores(i))
            {
              min_score=scores(i);
              min_index=i;
            }
          cout << "min_index = " << min_index << endl;
        }
      
      //The case where we need to go more steps
      if(min_index==steps-1)
        {
          cout << "Bracket Stepping Further --------------" << endl;
          //make gradient magnitude bigger
          scale=scale*pow(2.0,(double)(steps-3));
          //try to bracket starting with the second to last index, this 
          //way we know

          //that the first step will be downhill
          bracket(OF, position,searchDirection,alphas(steps-2),scale,
                  recursion_depth+1);
        }
      
      //The case where the first step was already too big
      //if(min is the first index)
      else if(min_index==0)
        {         
          cout << "Bracket Stepping Closer --------------" << endl;

          //make search direction vector smaller
          scale=scale*pow(2,-1.0*(steps-2));
          
          //try to bracket with smaller steps from x0
          bracket(OF,position,searchDirection,alphas(0),scale,
                  recursion_depth+1);
        }
      
      //The case where we have bracketed a min
      //if(min is not the first or last index) means we have bracketed the min
      else
        {
          //make current alphas and scores updated to reflect the tighter 
          //bracket
          alpha_left   = alphas(min_index-1);
          score_left   = scores(min_index-1);
          alpha_middle = alphas(min_index);
          score_middle = scores(min_index);
          alpha_right  = alphas(min_index+1);
          score_right  = scores(min_index+1);
        }
    }
  else
    {
      alpha_left   = -1.0;
      alpha_middle = -1.0;
      alpha_right  = -1.0;
    }
}

//This will make a properly bracketed set of alphas tighter around a minimum
void QMCMikesBracketingStepLengthSelector::quadratic_rebracketer(
  QMCObjectiveFunction * OF, Array1D<double> & position, Array1D<double> & searchDirection)
{
  //how many trial points will we try to rebracket
  int trials=4;

  //allocate a population to be evaulated
  Array1D < Array1D <double> > population;
  population.allocate(trials);
  for(int i=0;i<trials;i++)
    {
      population(i).allocate(position.dim1());
    }

  //allocate the alphas to be tried
  Array1D<double> alphas;
  alphas.allocate(trials);

  //allocate the scores to be tried
  Array1D<double> scores;
  scores.allocate(trials);

  //allocate the results from an objective function evaulation
  Array1D< QMCObjectiveFunctionResult > results;
  results.allocate(trials);

  //the one alpha which is the most promising will be called alpha_star
  double alpha_star;

  //determine minimum (alpha_star) of parabola with current alphas and scores
  double a=alpha_left;
  double b=alpha_middle;
  double c=alpha_right;
  double fa=score_left;
  double fb=score_middle;
  double fc=score_right;
  double bafbfc=(b-a)*(fb-fc);  
  double bcfbfa=(b-c)*(fb-fa);

  //actual calculate this minimum Numerical Recipes in C page 402 (10.2)
  alpha_star=b-0.5*((b-a)*bafbfc-(b-c)*bcfbfa)/(bafbfc-bcfbfa);

  //establish measures to help push the alpha_star away from current points
  double d_LM=alpha_middle-alpha_left;
  double d_MR=alpha_right-alpha_middle; 

  //how close we let alpha_star get to a currently calculated point
  double buffer_fraction=0.1; 

  //if alpha_star is too close to alpha_left push it away
  if(alpha_star<(alpha_left+buffer_fraction*d_LM))
    {
      alpha_star=alpha_left+buffer_fraction*d_LM;
    }

  //if alpha_star is too close to alpha_middle push it away
   if(alpha_star>(alpha_middle-buffer_fraction*d_LM))
    {
      alpha_star=alpha_middle-buffer_fraction*d_LM;
    }
   if(alpha_star<(alpha_middle+buffer_fraction*d_MR))
    {
      alpha_star=alpha_middle+buffer_fraction*d_MR;
    } 

  //if alpha_star is too close to alpha_right push it away
  if(alpha_star>(alpha_right-buffer_fraction*d_MR))
    {
      alpha_star=alpha_right-buffer_fraction*d_MR;
    }

  //put the alpha_star into the alpha vector and 
  //put the other three guess alphas in the alpha vector
  if(alpha_star<alpha_middle)
    {
      alphas(0) = (alpha_left+alpha_star)/2.0;
      alphas(1) = alpha_star;
      alphas(2) = (alpha_star+alpha_middle)/2.0;
      alphas(3) = (alpha_middle+alpha_right)/2.0;
    }
  else
    {
      alphas(0) = (alpha_left+alpha_middle)/2.0;
      alphas(1) = (alpha_middle+alpha_star)/2.0;
      alphas(2) = alpha_star;
      alphas(3) = (alpha_star+alpha_right)/2.0;
    }

  //put the guess vectors in the population
  for(int i=0;i<population.dim1();i++)
    {
      for(int j=0;j<population(i).dim1();j++)
        {
          population(i)(j)=position(j)+alphas(i)*
            searchDirection(j);
        }
    }

  //evaluate the population
  results = OF->evaluate(population);

  //put the scores into the scores vector
  for(int i=0;i<population.dim1();i++)
    {
      scores(i)=results(i).getScore();
    }

  //find the minimum in {population,alpha_middle} 
  //(global min of current points)
  int min_index=-1; //if at the end the min index is -1 we 
                    //know the old middle value is still the best score
  double min_score=score_middle;
  for(int i=0;i<population.dim1();i++)
    {
      if(min_score>scores(i))
        {
          min_score=scores(i);
          min_index=i;
        }
    }

  //make current alphas and scores updated to reflect the tighter bracket
  if(min_index==-1)   
    {   
      //the old min is still the min so keep the "middle" where it is
      if(alpha_star<alpha_middle)
        {  
          // {left,0,star(1),2,middle,3,right}
          alpha_left  = alphas(2);
          score_left  = scores(2);
          alpha_right = alphas(3);
          score_right = scores(3);
        }
      else
        { 
          // {left,0,middle,1,star(2),3,right}
          alpha_left  = alphas(0);
          score_left  = scores(0);
          alpha_right = alphas(1);
          score_right = scores(1);
        }
    }
  else
    {
      if(alpha_star<alpha_middle)
        { 
          // {left,0,star(1),2,middle,3,right}
          if(min_index==0)
            {  
              alpha_middle = alphas(0);
              score_middle = scores(0);
              alpha_right  = alphas(1);
              score_right  = scores(1);
            }
          else if(min_index==1)
            { 
              alpha_left   = alphas(0);
              score_left   = scores(0);
              alpha_middle = alphas(1);
              score_middle = scores(1);
              alpha_right  = alphas(2);
              score_right  = scores(2);
            }
          else if(min_index==2)
            { 
              alpha_right  = alpha_middle;
              score_right  = score_middle;
              alpha_left   = alphas(1);
              score_left   = scores(1);
              alpha_middle = alphas(2);
              score_middle = scores(2);
            }
          else
            { 
              alpha_left   = alpha_middle;
              score_left   = score_middle;
              alpha_middle = alphas(3);
              score_middle = scores(3);
            }
        }
      else
        { 
          // {left,0,middle,1,star(2),3,right}
          if(min_index==0)
            {  
              alpha_right  = alpha_middle;
              score_right  = score_middle;
              alpha_middle = alphas(0);
              score_middle = scores(0);
            }
          else if(min_index==1)
            { 
              alpha_left   = alpha_middle;
              score_left   = score_middle;
              alpha_middle = alphas(1);
              score_middle = scores(1);
              alpha_right  = alphas(2);
              score_right  = scores(2);
            }
          else if(min_index==2)
            { 
              alpha_left   = alphas(1);
              score_left   = scores(1);
              alpha_middle = alphas(2);
              score_middle = scores(2);
              alpha_right  = alphas(3);
              score_right  = scores(3);
            }
          else
            { 
              alpha_left   = alphas(2);
              score_left   = scores(2);
              alpha_middle = alphas(3);
              score_middle = scores(3);
            }  
        }
    }
}















