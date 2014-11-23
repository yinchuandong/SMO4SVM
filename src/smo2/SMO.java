package smo2;

/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    SMO.java
 *    Copyright (C) 1999 Eibe Frank
 *
 */

//package weka.classifiers.functions;
import weka.classifiers.functions.supportVector.*;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.Logistic;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
import weka.filters.unsupervised.attribute.Normalize;
import weka.filters.unsupervised.attribute.Standardize;
import weka.filters.Filter;

import java.util.*;
import java.lang.Math;
import java.io.*;

import weka.core.*;
@SuppressWarnings({"unchecked","serial","cast"})
/**
 * Implements John C. Platt's sequential minimal optimization
 * algorithm for training a support vector classifier using polynomial
 * or RBF kernels. 
 *
 * This implementation globally replaces all missing values and
 * transforms nominal attributes into binary ones. It also
 * normalizes all attributes by default. (Note that the coefficients
 * in the output are based on the normalized/standardized data, not the
 * original data.)
 *
 * Multi-class problems are solved using pairwise classification.
 *
 * To obtain proper probability estimates, use the option that fits
 * logistic regression models to the outputs of the support vector
 * machine. In the multi-class case the predicted probabilities
 * will be coupled using Hastie and Tibshirani's pairwise coupling
 * method.
 *
 * Note: for improved speed standardization should be turned off when
 * operating on SparseInstances.<p>
 *
 * For more information on the SMO algorithm, see<p>
 *
 * J. Platt (1998). <i>Fast Training of Support Vector
 * Machines using Sequential Minimal Optimization</i>. Advances in Kernel
 * Methods - Support Vector Learning, B. Schoelkopf, C. Burges, and
 * A. Smola, eds., MIT Press. <p>
 *
 * S.S. Keerthi, S.K. Shevade, C. Bhattacharyya, K.R.K. Murthy, 
 * <i>Improvements to Platt's SMO Algorithm for SVM Classifier Design</i>. 
 * Neural Computation, 13(3), pp 637-649, 2001. <p>
 *
 * Valid options are:<p>
 *
 * -C num <br>
 * The complexity constant C. (default 1)<p>
 *
 * -E num <br>
 * The exponent for the polynomial kernel. (default 1)<p>
 *
 * -G num <br>
 * Gamma for the RBF kernel. (default 0.01)<p>
 *
 * -N <0|1|2> <br>
 * Whether to 0=normalize/1=standardize/2=neither. (default 0=normalize)<p>
 *
 * -F <br>
 * Feature-space normalization (only for non-linear polynomial kernels). <p>
 *
 * -O <br>
 * Use lower-order terms (only for non-linear polynomial kernels). <p>
 *
 * -R <br>
 * Use the RBF kernel. (default poly)<p>
 *
 * -A num <br>
 * Sets the size of the kernel cache. Should be a prime number. 
 * (default 250007, use 0 for full cache) <p>
 *
 * -L num <br>
 * Sets the tolerance parameter. (default 1.0e-3)<p>
 *
 * -P num <br>
 * Sets the epsilon for round-off error. (default 1.0e-12)<p>
 *
 * -M <br>
 * Fit logistic models to SVM outputs.<p>
 *
 * -V num <br>
 * Number of folds for cross-validation used to generate data
 * for logistic models. (default -1, use training data)
 *
 * -W num <br>
 * Random number seed for cross-validation. (default 1)
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Shane Legg (shane@intelligenesis.net) (sparse vector code)
 * @author Stuart Inglis (stuart@reeltwo.com) (sparse vector code)
 * @version $Revision: 1.53.2.2 $ */
public class SMO extends Classifier implements WeightedInstancesHandler {

  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "Implements John Platt's sequential minimal optimization "
      + "algorithm for training a support vector classifier.\n\n"
      + "This implementation globally replaces all missing values and "
      + "transforms nominal attributes into binary ones. It also "
      + "normalizes all attributes by default. (In that case the coefficients "
      + "in the output are based on the normalized data, not the "
      + "original data --- this is important for interpreting the classifier.)\n\n"
      + "Multi-class problems are solved using pairwise classification.\n\n"
      + "To obtain proper probability estimates, use the option that fits "
      + "logistic regression models to the outputs of the support vector "
      + "machine. In the multi-class case the predicted probabilities "
      + "are coupled using Hastie and Tibshirani's pairwise coupling "
      + "method.\n\n"
      + "Note: for improved speed normalization should be turned off when "
      + "operating on SparseInstances.\n\n"
      + "For more information on the SMO algorithm, see\n\n"
      + "J. Platt (1998). \"Fast Training of Support Vector "
      + "Machines using Sequential Minimal Optimization\". Advances in Kernel "
      + "Methods - Support Vector Learning, B. Schoelkopf, C. Burges, and "
      + "A. Smola, eds., MIT Press. \n\n"
      + "S.S. Keerthi, S.K. Shevade, C. Bhattacharyya, K.R.K. Murthy,  "
      + "\"Improvements to Platt's SMO Algorithm for SVM Classifier Design\".  "
      + "Neural Computation, 13(3), pp 637-649, 2001.";
  }

  /**
   * Class for building a binary support vector machine.
   */
  protected class BinarymySMO implements Serializable {
    
    /** The Lagrange multipliers. */
    protected double[] m_alpha;

    /** The thresholds. */
    protected double m_b, m_bLow, m_bUp;

    /** The indices for m_bLow and m_bUp */
    protected int m_iLow, m_iUp;

    /** The training data. */
    protected Instances m_data;

    /** Weight vector for linear machine. */
    protected double[] m_weights;

    /** Variables to hold weight vector in sparse form.
	(To reduce storage requirements.) */
    protected double[] m_sparseWeights;
    protected int[] m_sparseIndices;

    /** Kernel to use **/
    protected Kernel m_kernel;

    /** The transformed class values. */
    protected double[] m_class;

    /** The current set of errors for all non-bound examples. */
    protected double[] m_errors;

    /** The five different sets used by the algorithm. */
    protected SMOset m_I0; // {i: 0 < m_alpha[i] < C}
    protected SMOset m_I1; // {i: m_class[i] = 1, m_alpha[i] = 0}
    protected SMOset m_I2; // {i: m_class[i] = -1, m_alpha[i] =C}
    protected SMOset m_I3; // {i: m_class[i] = 1, m_alpha[i] = C}
    protected SMOset m_I4; // {i: m_class[i] = -1, m_alpha[i] = 0}

    /** The set of support vectors */
    protected SMOset m_supportVectors; // {i: 0 < m_alpha[i]}

    /** Stores logistic regression model for probability estimate */
    protected Logistic m_logistic = null;

    /** Stores the weight of the training instances */
    protected double m_sumOfWeights = 0;

    /**
     * Fits logistic regression model to SVM outputs analogue
     * to John Platt's method.  
     *
     * @param insts the set of training instances
     * @param cl1 the first class' index
     * @param cl2 the second class' index
     * @exception Exception if the sigmoid can't be fit successfully
     */
    protected void fitLogistic(Instances insts, int cl1, int cl2,
			     int numFolds, Random random) 
      throws Exception {

      // Create header of instances object
      FastVector atts = new FastVector(2);
      atts.addElement(new Attribute("pred"));
      FastVector attVals = new FastVector(2);
      attVals.addElement(insts.classAttribute().value(cl1));
      attVals.addElement(insts.classAttribute().value(cl2));
      atts.addElement(new Attribute("class", attVals));
      Instances data = new Instances("data", atts, insts.numInstances());
      data.setClassIndex(1);

      // Collect data for fitting the logistic model
      if (numFolds <= 0) {

	// Use training data
	for (int j = 0; j < insts.numInstances(); j++) {
	  Instance inst = insts.instance(j);
	  double[] vals = new double[2];
	  vals[0] = mySVMOutput(-1, inst);
	  if (inst.classValue() == cl2) {
	    vals[1] = 1;
	  }
	  data.add(new Instance(inst.weight(), vals));
	}
      } else {

	// Check whether number of folds too large
	if (numFolds > insts.numInstances()) {
	  numFolds = insts.numInstances();
	}

	// Make copy of instances because we will shuffle them around
	insts = new Instances(insts);
	
	// Perform three-fold cross-validation to collect
	// unbiased predictions
	insts.randomize(random);
	insts.stratify(numFolds);
	for (int i = 0; i < numFolds; i++) {
		//这一句版本有问题
//	  Instances train = insts.trainCV(numFolds, i, random);
	Instances train = insts.trainCV(numFolds, i);
	  SerializedObject so = new SerializedObject(this);
	  BinarymySMO mySMO = (BinarymySMO)so.getObject();
	  mySMO.buildClassifier(train, cl1, cl2, false, -1, -1);
	  Instances test = insts.testCV(numFolds, i);
	  for (int j = 0; j < test.numInstances(); j++) {
	    double[] vals = new double[2];
	    vals[0] = mySMO.mySVMOutput(-1, test.instance(j));
	    if (test.instance(j).classValue() == cl2) {
	      vals[1] = 1;
	    }
	    data.add(new Instance(test.instance(j).weight(), vals));
	  }
	}
      }

      // Build logistic regression model
      m_logistic = new Logistic();
      m_logistic.buildClassifier(data);
    }

    /**
     * Method for building the binary classifier.
     *
     * @param insts the set of training instances
     * @param cl1 the first class' index
     * @param cl2 the second class' index
     * @param fitLogistic true if logistic model is to be fit
     * @param numFolds number of folds for internal cross-validation
     * @param random random number generator for cross-validation
     * @exception Exception if the classifier can't be built successfully
     */
    protected void buildClassifier(Instances insts, int cl1, int cl2,
				 boolean fitLogistic, int numFolds,
				 int randomSeed) throws Exception {
      
      // Initialize some variables
      m_bUp = -1; m_bLow = 1; m_b = 0; 
      m_alpha = null; m_data = null; m_weights = null; m_errors = null;
      m_logistic = null; m_I0 = null; m_I1 = null; m_I2 = null;
      m_I3 = null; m_I4 = null;	m_sparseWeights = null; m_sparseIndices = null;

      // Store the sum of weights
      m_sumOfWeights = insts.sumOfWeights();
      
      // Set class values
      m_class = new double[insts.numInstances()];
      m_iUp = -1; m_iLow = -1;
      for (int i = 0; i < m_class.length; i++) {
	if ((int) insts.instance(i).classValue() == cl1) {
	  m_class[i] = -1; m_iLow = i;
	} else if ((int) insts.instance(i).classValue() == cl2) {
	  m_class[i] = 1; m_iUp = i;
	} else {
	  throw new Exception ("This should never happen!");
	}
      }

      // Check whether one or both classes are missing
      if ((m_iUp == -1) || (m_iLow == -1)) {
	if (m_iUp != -1) {
	  m_b = -1;
	} else if (m_iLow != -1) {
	  m_b = 1;
	} else {
	  m_class = null;
	  return;
	}
	if (!m_useRBF && m_exponent == 1.0) {
	  m_sparseWeights = new double[0];
	  m_sparseIndices = new int[0];
	  m_class = null;
	} else {
	  m_supportVectors = new SMOset(0);
	  m_alpha = new double[0];
	  m_class = new double[0];
	}

	// Fit sigmoid if requested
	if (fitLogistic) {
	  fitLogistic(insts, cl1, cl2, numFolds, new Random(randomSeed));
	}
	return;
      }
      
      // Set the reference to the data
      m_data = insts;

      // If machine is linear, reserve space for weights
      if (!m_useRBF && m_exponent == 1.0) {
	m_weights = new double[m_data.numAttributes()];
      } else {
	m_weights = null;
      }
      
      // Initialize alpha array to zero
      m_alpha = new double[m_data.numInstances()];
      
      // Initialize sets
      m_supportVectors = new SMOset(m_data.numInstances());
      m_I0 = new SMOset(m_data.numInstances());
      m_I1 = new SMOset(m_data.numInstances());
      m_I2 = new SMOset(m_data.numInstances());
      m_I3 = new SMOset(m_data.numInstances());
      m_I4 = new SMOset(m_data.numInstances());

      // Clean out some instance variables
      m_sparseWeights = null;
      m_sparseIndices = null;
      
      // Initialize error cache
      m_errors = new double[m_data.numInstances()];
      m_errors[m_iLow] = 1; m_errors[m_iUp] = -1;
     
      // Initialize kernel
      if(m_useRBF) {
	m_kernel = new RBFKernel(m_data, m_cacheSize, m_gamma);
      } else {
	if (m_featureSpaceNormalization) {
	  m_kernel = new NormalizedPolyKernel(m_data, m_cacheSize, m_exponent, 
					      m_lowerOrder);
	} else {
	  m_kernel = new PolyKernel(m_data, m_cacheSize, m_exponent, m_lowerOrder);
	}
      }
      
      // Build up I1 and I4
      for (int i = 0; i < m_class.length; i++ ) {
	if (m_class[i] == 1) {
	  m_I1.insert(i);
	} else {
	  m_I4.insert(i);
	}
      }
      
      // Loop to find all the support vectors
      int numChanged = 0;
      boolean examineAll = true;
      while ((numChanged > 0) || examineAll) {
	numChanged = 0;
	if (examineAll) {
	  for (int i = 0; i < m_alpha.length; i++) {
	    if (examineExample(i)) {
	      numChanged++;
	    }
	  }
	} else {
	  
	  // This code implements Modification 1 from Keerthi et al.'s paper
	  for (int i = 0; i < m_alpha.length; i++) {
	    if ((m_alpha[i] > 0) &&  
		(m_alpha[i] < m_C * m_data.instance(i).weight())) {
	      if (examineExample(i)) {
		numChanged++;
	      }
	      
	      // Is optimality on unbound vectors obtained?
	      if (m_bUp > m_bLow - 2 * m_tol) {
		numChanged = 0;
		break;
	      }
	    }
	  }
	  
	  //This is the code for Modification 2 from Keerthi et al.'s paper
	  /*boolean innerLoopSuccess = true; 
	    numChanged = 0;
	    while ((m_bUp < m_bLow - 2 * m_tol) && (innerLoopSuccess == true)) {
	    innerLoopSuccess = takeStep(m_iUp, m_iLow, m_errors[m_iLow]);
	    }*/
	}
	
	if (examineAll) {
	  examineAll = false;
	} else if (numChanged == 0) {
	  examineAll = true;
	}
      }
      
      // Set threshold
      m_b = (m_bLow + m_bUp) / 2.0;
      
      // Save memory
      m_kernel.clean(); 
      
      m_errors = null;
      m_I0 = m_I1 = m_I2 = m_I3 = m_I4 = null;
      
      // If machine is linear, delete training data
      // and store weight vector in sparse format
      if (!m_useRBF && m_exponent == 1.0) {
	
	// We don't need to store the set of support vectors
	m_supportVectors = null;

	// We don't need to store the class values either
	m_class = null;
	
	// Clean out training data
	if (!m_checksTurnedOff) {
	  m_data = new Instances(m_data, 0);
	} else {
	  m_data = null;
	}
	
	// Convert weight vector
	double[] sparseWeights = new double[m_weights.length];
	int[] sparseIndices = new int[m_weights.length];
	int counter = 0;
	for (int i = 0; i < m_weights.length; i++) {
	  if (m_weights[i] != 0.0) {
	    sparseWeights[counter] = m_weights[i];
	    sparseIndices[counter] = i;
	    counter++;
	  }
	}
	m_sparseWeights = new double[counter];
	m_sparseIndices = new int[counter];
	System.arraycopy(sparseWeights, 0, m_sparseWeights, 0, counter);
	System.arraycopy(sparseIndices, 0, m_sparseIndices, 0, counter);
	
	// Clean out weight vector
	m_weights = null;
	
	// We don't need the alphas in the linear case
	m_alpha = null;
      }
      
      // Fit sigmoid if requested
      if (fitLogistic) {
	fitLogistic(insts, cl1, cl2, numFolds, new Random(randomSeed));
      }

    }
    
    /**
     * Computes SVM output for given instance.
     *
     * @param index the instance for which output is to be computed
     * @param inst the instance 
     * @return the output of the SVM for the given instance
     */
    protected double mySVMOutput(int index, Instance inst) throws Exception {
      
      double result = 0;
      
      // Is the machine linear?
      if (!m_useRBF && m_exponent == 1.0) {
	
	// Is weight vector stored in sparse format?
	if (m_sparseWeights == null) {
	  int n1 = inst.numValues(); 
	  for (int p = 0; p < n1; p++) {
	    if (inst.index(p) != m_classIndex) {
	      result += m_weights[inst.index(p)] * inst.valueSparse(p);
	    }
	  }
	} else {
	  int n1 = inst.numValues(); int n2 = m_sparseWeights.length;
	  for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
	    int ind1 = inst.index(p1); 
	    int ind2 = m_sparseIndices[p2];
	    if (ind1 == ind2) {
	      if (ind1 != m_classIndex) {
		result += inst.valueSparse(p1) * m_sparseWeights[p2];
	      }
	      p1++; p2++;
	    } else if (ind1 > ind2) {
	      p2++;
	    } else { 
	      p1++;
	    }
	  }
	}
      } else {
	for (int i = m_supportVectors.getNext(-1); i != -1; 
	     i = m_supportVectors.getNext(i)) {
	  result += m_class[i] * m_alpha[i] * m_kernel.eval(index, i, inst);
	}
      }
      result -= m_b;
      
      return result;
    }

    /**
     * Prints out the classifier.
     *
     * @return a description of the classifier as a string
     */
    public String toString() {

      StringBuffer text = new StringBuffer();
      int printed = 0;

      if ((m_alpha == null) && (m_sparseWeights == null)) {
	return "BinarymySMO: No model built yet.\n";
      }
      try {
	text.append("BinarymySMO\n\n");

	// If machine linear, print weight vector
	if (!m_useRBF && m_exponent == 1.0) {
	  text.append("Machine linear: showing attribute weights, ");
	  text.append("not support vectors.\n\n");

	  // We can assume that the weight vector is stored in sparse
	  // format because the classifier has been built
	  for (int i = 0; i < m_sparseWeights.length; i++) {
	    if (m_sparseIndices[i] != (int)m_classIndex) {
	      if (printed > 0) {
		text.append(" + ");
	      } else {
		text.append("   ");
	      }
	      text.append(Utils.doubleToString(m_sparseWeights[i], 12, 4) +
			  " * ");
	      if (m_filterType == FILTER_STANDARDIZE) {
		text.append("(standardized) ");
	      } else if (m_filterType == FILTER_NORMALIZE) {
		text.append("(normalized) ");
	      }
	      if (!m_checksTurnedOff) {
		text.append(m_data.attribute(m_sparseIndices[i]).name()+"\n");
	      } else {
		text.append("attribute with index " + 
			    m_sparseIndices[i] +"\n");
	      }
	      printed++;
	    }
	  }
	} else {
	  for (int i = 0; i < m_alpha.length; i++) {
	    if (m_supportVectors.contains(i)) {
	      double val = m_alpha[i];
	      if (m_class[i] == 1) {
		if (printed > 0) {
		  text.append(" + ");
		}
	      } else {
		text.append(" - ");
	      }
	      text.append(Utils.doubleToString(val, 12, 4) 
			  + " * <");
	      for (int j = 0; j < m_data.numAttributes(); j++) {
		if (j != m_data.classIndex()) {
		  text.append(m_data.instance(i).toString(j));
		}
		if (j != m_data.numAttributes() - 1) {
		  text.append(" ");
		}
	      }
	      text.append("> * X]\n");
	      printed++;
	    }
	  }
	}
	if (m_b > 0) {
	  text.append(" - " + Utils.doubleToString(m_b, 12, 4));
	} else {
	  text.append(" + " + Utils.doubleToString(-m_b, 12, 4));
	}

	if (m_useRBF || m_exponent != 1.0) {
	  text.append("\n\nNumber of support vectors: " + 
		      m_supportVectors.numElements());
	}
	int numEval = 0;
	int numCacheHits = -1;
	if(m_kernel != null)
	{
	  numEval = m_kernel.numEvals();
	  //这一句有问题，没有改函数
//	  numCacheHits = m_kernel.numCacheHits();
	}
	text.append("\n\nNumber of kernel evaluations: " + numEval);
	if (numCacheHits >= 0 && numEval > 0)
	{
		double hitRatio = 1 - numEval*1.0/(numCacheHits+numEval);
		text.append(" (" + Utils.doubleToString(hitRatio*100, 7, 3).trim() + "% cached)");
	}

      } catch (Exception e) {
	e.printStackTrace();

	return "Can't print BinarymySMO classifier.";
      }
    
      return text.toString();
    }

    /**
     * Examines instance.
     *
     * @param i2 index of instance to examine
     * @return true if examination was successfull
     * @exception Exception if something goes wrong
     */
    protected boolean examineExample(int i2) throws Exception {
    
      double y2, alph2, F2;
      int i1 = -1;
    
      y2 = m_class[i2];
      alph2 = m_alpha[i2];
      if (m_I0.contains(i2)) {
	F2 = m_errors[i2];
      } else {
	F2 = mySVMOutput(i2, m_data.instance(i2)) + m_b - y2;
	m_errors[i2] = F2;
      
	// Update thresholds
	if ((m_I1.contains(i2) || m_I2.contains(i2)) && (F2 < m_bUp)) {
	  m_bUp = F2; m_iUp = i2;
	} else if ((m_I3.contains(i2) || m_I4.contains(i2)) && (F2 > m_bLow)) {
	  m_bLow = F2; m_iLow = i2;
	}
      }

      // Check optimality using current bLow and bUp and, if
      // violated, find an index i1 to do joint optimization
      // with i2...
      boolean optimal = true;
      if (m_I0.contains(i2) || m_I1.contains(i2) || m_I2.contains(i2)) {
	if (m_bLow - F2 > 2 * m_tol) {
	  optimal = false; i1 = m_iLow;
	}
      }
      if (m_I0.contains(i2) || m_I3.contains(i2) || m_I4.contains(i2)) {
	if (F2 - m_bUp > 2 * m_tol) {
	  optimal = false; i1 = m_iUp;
	}
      }
      if (optimal) {
	return false;
      }

      // For i2 unbound choose the better i1...
      if (m_I0.contains(i2)) {
	if (m_bLow - F2 > F2 - m_bUp) {
	  i1 = m_iLow;
	} else {
	  i1 = m_iUp;
	}
      }
      if (i1 == -1) {
	throw new Exception("This should never happen!");
      }
      return takeStep(i1, i2, F2);
    }

    /**
     * Method solving for the Lagrange multipliers for
     * two instances.
     *
     * @param i1 index of the first instance
     * @param i2 index of the second instance
     * @return true if multipliers could be found
     * @exception Exception if something goes wrong
     */
    protected boolean takeStep(int i1, int i2, double F2) throws Exception {

      double alph1, alph2, y1, y2, F1, s, L, H, k11, k12, k22, eta,
	a1, a2, f1, f2, v1, v2, Lobj, Hobj, b1, b2, bOld;
      double C1 = m_C * m_data.instance(i1).weight();
      double C2 = m_C * m_data.instance(i2).weight();

      // Don't do anything if the two instances are the same
      if (i1 == i2) {
	return false;
      }

      // Initialize variables
      alph1 = m_alpha[i1]; alph2 = m_alpha[i2];
      y1 = m_class[i1]; y2 = m_class[i2];
      F1 = m_errors[i1];
      s = y1 * y2;

      // Find the constraints on a2
      if (y1 != y2) {
	L = Math.max(0, alph2 - alph1); 
	H = Math.min(C2, C1 + alph2 - alph1);
      } else {
	L = Math.max(0, alph1 + alph2 - C1);
	H = Math.min(C2, alph1 + alph2);
      }
      if (L >= H) {
	return false;
      }

      // Compute second derivative of objective function
      k11 = m_kernel.eval(i1, i1, m_data.instance(i1));
      k12 = m_kernel.eval(i1, i2, m_data.instance(i1));
      k22 = m_kernel.eval(i2, i2, m_data.instance(i2));
      eta = 2 * k12 - k11 - k22;

      // Check if second derivative is negative
      if (eta < 0) {

	// Compute unconstrained maximum
	a2 = alph2 - y2 * (F1 - F2) / eta;

	// Compute constrained maximum
	if (a2 < L) {
	  a2 = L;
	} else if (a2 > H) {
	  a2 = H;
	}
      } else {

	// Look at endpoints of diagonal
	f1 = mySVMOutput(i1, m_data.instance(i1));
	f2 = mySVMOutput(i2, m_data.instance(i2));
	v1 = f1 + m_b - y1 * alph1 * k11 - y2 * alph2 * k12; 
	v2 = f2 + m_b - y1 * alph1 * k12 - y2 * alph2 * k22; 
	double gamma = alph1 + s * alph2;
	Lobj = (gamma - s * L) + L - 0.5 * k11 * (gamma - s * L) * (gamma - s * L) - 
	  0.5 * k22 * L * L - s * k12 * (gamma - s * L) * L - 
	  y1 * (gamma - s * L) * v1 - y2 * L * v2;
	Hobj = (gamma - s * H) + H - 0.5 * k11 * (gamma - s * H) * (gamma - s * H) - 
	  0.5 * k22 * H * H - s * k12 * (gamma - s * H) * H - 
	  y1 * (gamma - s * H) * v1 - y2 * H * v2;
	if (Lobj > Hobj + m_eps) {
	  a2 = L;
	} else if (Lobj < Hobj - m_eps) {
	  a2 = H;
	} else {
	  a2 = alph2;
	}
      }
      if (Math.abs(a2 - alph2) < m_eps * (a2 + alph2 + m_eps)) {
	return false;
      }
      
      // To prevent precision problems
      if (a2 > C2 - m_Del * C2) {
	a2 = C2;
      } else if (a2 <= m_Del * C2) {
	a2 = 0;
      }
      
      // Recompute a1
      a1 = alph1 + s * (alph2 - a2);
      
      // To prevent precision problems
      if (a1 > C1 - m_Del * C1) {
	a1 = C1;
      } else if (a1 <= m_Del * C1) {
	a1 = 0;
      }
      
      // Update sets
      if (a1 > 0) {
	m_supportVectors.insert(i1);
      } else {
	m_supportVectors.delete(i1);
      }
      if ((a1 > 0) && (a1 < C1)) {
	m_I0.insert(i1);
      } else {
	m_I0.delete(i1);
      }
      if ((y1 == 1) && (a1 == 0)) {
	m_I1.insert(i1);
      } else {
	m_I1.delete(i1);
      }
      if ((y1 == -1) && (a1 == C1)) {
	m_I2.insert(i1);
      } else {
	m_I2.delete(i1);
      }
      if ((y1 == 1) && (a1 == C1)) {
	m_I3.insert(i1);
      } else {
	m_I3.delete(i1);
      }
      if ((y1 == -1) && (a1 == 0)) {
	m_I4.insert(i1);
      } else {
	m_I4.delete(i1);
      }
      if (a2 > 0) {
	m_supportVectors.insert(i2);
      } else {
	m_supportVectors.delete(i2);
      }
      if ((a2 > 0) && (a2 < C2)) {
	m_I0.insert(i2);
      } else {
	m_I0.delete(i2);
      }
      if ((y2 == 1) && (a2 == 0)) {
	m_I1.insert(i2);
      } else {
	m_I1.delete(i2);
      }
      if ((y2 == -1) && (a2 == C2)) {
	m_I2.insert(i2);
      } else {
	m_I2.delete(i2);
      }
      if ((y2 == 1) && (a2 == C2)) {
	m_I3.insert(i2);
      } else {
	m_I3.delete(i2);
      }
      if ((y2 == -1) && (a2 == 0)) {
	m_I4.insert(i2);
      } else {
	m_I4.delete(i2);
      }
      
      // Update weight vector to reflect change a1 and a2, if linear SVM
      if (!m_useRBF && m_exponent == 1.0) {
	Instance inst1 = m_data.instance(i1);
	for (int p1 = 0; p1 < inst1.numValues(); p1++) {
	  if (inst1.index(p1) != m_data.classIndex()) {
	    m_weights[inst1.index(p1)] += 
	      y1 * (a1 - alph1) * inst1.valueSparse(p1);
	  }
	}
	Instance inst2 = m_data.instance(i2);
	for (int p2 = 0; p2 < inst2.numValues(); p2++) {
	  if (inst2.index(p2) != m_data.classIndex()) {
	    m_weights[inst2.index(p2)] += 
	      y2 * (a2 - alph2) * inst2.valueSparse(p2);
	  }
	}
      }
      
      // Update error cache using new Lagrange multipliers
      for (int j = m_I0.getNext(-1); j != -1; j = m_I0.getNext(j)) {
	if ((j != i1) && (j != i2)) {
	  m_errors[j] += 
	    y1 * (a1 - alph1) * m_kernel.eval(i1, j, m_data.instance(i1)) + 
	    y2 * (a2 - alph2) * m_kernel.eval(i2, j, m_data.instance(i2));
	}
      }
      
      // Update error cache for i1 and i2
      m_errors[i1] += y1 * (a1 - alph1) * k11 + y2 * (a2 - alph2) * k12;
      m_errors[i2] += y1 * (a1 - alph1) * k12 + y2 * (a2 - alph2) * k22;
      
      // Update array with Lagrange multipliers
      m_alpha[i1] = a1;
      m_alpha[i2] = a2;
      
      // Update thresholds
      m_bLow = -Double.MAX_VALUE; m_bUp = Double.MAX_VALUE;
      m_iLow = -1; m_iUp = -1;
      for (int j = m_I0.getNext(-1); j != -1; j = m_I0.getNext(j)) {
	if (m_errors[j] < m_bUp) {
	  m_bUp = m_errors[j]; m_iUp = j;
	}
	if (m_errors[j] > m_bLow) {
	  m_bLow = m_errors[j]; m_iLow = j;
	}
      }
      if (!m_I0.contains(i1)) {
	if (m_I3.contains(i1) || m_I4.contains(i1)) {
	  if (m_errors[i1] > m_bLow) {
	    m_bLow = m_errors[i1]; m_iLow = i1;
	  } 
	} else {
	  if (m_errors[i1] < m_bUp) {
	    m_bUp = m_errors[i1]; m_iUp = i1;
	  }
	}
      }
      if (!m_I0.contains(i2)) {
	if (m_I3.contains(i2) || m_I4.contains(i2)) {
	  if (m_errors[i2] > m_bLow) {
	    m_bLow = m_errors[i2]; m_iLow = i2;
	  }
	} else {
	  if (m_errors[i2] < m_bUp) {
	    m_bUp = m_errors[i2]; m_iUp = i2;
	  }
	}
      }
      if ((m_iLow == -1) || (m_iUp == -1)) {
	throw new Exception("This should never happen!");
      }

      // Made some progress.
      return true;
    }
  
    /**
     * Quick and dirty check whether the quadratic programming problem is solved.
     */
    protected void checkClassifier() throws Exception {

      double sum = 0;
      for (int i = 0; i < m_alpha.length; i++) {
	if (m_alpha[i] > 0) {
	  sum += m_class[i] * m_alpha[i];
	}
      }
      System.err.println("Sum of y(i) * alpha(i): " + sum);

      for (int i = 0; i < m_alpha.length; i++) {
	double output = mySVMOutput(i, m_data.instance(i));
	if (Utils.eq(m_alpha[i], 0)) {
	  if (Utils.sm(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 1 violated: " + m_class[i] * output);
	  }
	} 
	if (Utils.gr(m_alpha[i], 0) && 
	    Utils.sm(m_alpha[i], m_C * m_data.instance(i).weight())) {
	  if (!Utils.eq(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 2 violated: " + m_class[i] * output);
	  }
	} 
	if (Utils.eq(m_alpha[i], m_C * m_data.instance(i).weight())) {
	  if (Utils.gr(m_class[i] * output, 1)) {
	    System.err.println("KKT condition 3 violated: " + m_class[i] * output);
	  }
	} 
      }
    }  
  }

  /** The filter to apply to the training data */
  public static final int FILTER_NORMALIZE = 0;
  public static final int FILTER_STANDARDIZE = 1;
  public static final int FILTER_NONE = 2;
  public static final Tag [] TAGS_FILTER = {
    new Tag(FILTER_NORMALIZE, "Normalize training data"),
    new Tag(FILTER_STANDARDIZE, "Standardize training data"),
    new Tag(FILTER_NONE, "No normalization/standardization"),
  };

  /** The binary classifier(s) */
  protected BinarymySMO[][] m_classifiers = null;

  /** The exponent for the polynomial kernel. */
  protected double m_exponent = 1.0;
 
  /** Use lower-order terms? */
  protected boolean m_lowerOrder = false;
  
  /** Gamma for the RBF kernel. */
  protected double m_gamma = 0.01;
  
  /** The complexity parameter. */
  protected double m_C = 1.0;
  
  /** Epsilon for rounding. */
  protected double m_eps = 1.0e-12;
  
  /** Tolerance for accuracy of result. */
  protected double m_tol = 1.0e-3;

  /** Whether to normalize/standardize/neither */
  protected int m_filterType = FILTER_NORMALIZE;
  
  /** Feature-space normalization? */
  protected boolean m_featureSpaceNormalization = false;
  
  /** Use RBF kernel? (default: poly) */
  protected boolean m_useRBF = false;
  
  /** The size of the cache (a prime number) */
  protected int m_cacheSize = 250007;

  /** The filter used to make attributes numeric. */
  protected NominalToBinary m_NominalToBinary;

  /** The filter used to standardize/normalize all values. */
  protected Filter m_Filter = null;

  /** The filter used to get rid of missing values. */
  protected ReplaceMissingValues m_Missing;

  /** Only numeric attributes in the dataset? */
  protected boolean m_onlyNumeric;

  /** The class index from the training data */
  protected int m_classIndex = -1;

  /** The class attribute */
  protected Attribute m_classAttribute;

  /** Turn off all checks and conversions? Turning them off assumes
      that data is purely numeric, doesn't contain any missing values,
      and has a nominal class. Turning them off also means that
      no header information will be stored if the machine is linear. 
      Finally, it also assumes that no instance has a weight equal to 0.*/
  protected boolean m_checksTurnedOff;

  /** Precision constant for updating sets */
  protected static double m_Del = 1000 * Double.MIN_VALUE;

  /** Whether logistic models are to be fit */
  protected boolean m_fitLogisticModels = false;

  /** The number of folds for the internal cross-validation */
  protected int m_numFolds = -1;

  /** The random number seed  */
  protected int m_randomSeed = 1;

  /**
   * Turns off checks for missing values, etc. Use with caution.
   */
  public void turnChecksOff() {

    m_checksTurnedOff = true;
  }

  /**
   * Turns on checks for missing values, etc.
   */
  public void turnChecksOn() {

    m_checksTurnedOff = false;
  }

  /**
   * Method for building the classifier. Implements a one-against-one
   * wrapper for multi-class problems.
   *
   * @param insts the set of training instances
   * @exception Exception if the classifier can't be built successfully
   */
  public void buildClassifier(Instances insts) throws Exception {

    if (!m_checksTurnedOff) {
      if (insts.checkForStringAttributes()) {
	throw new UnsupportedAttributeTypeException("Cannot handle string attributes!");
      }
      if (insts.classAttribute().isNumeric()) {
	throw new UnsupportedClassTypeException("mySMO can't handle a numeric class! Use"
						+ "SMOreg for performing regression.");
      }
      insts = new Instances(insts);
      insts.deleteWithMissingClass();
      if (insts.numInstances() == 0) {
	throw new Exception("No training instances without a missing class!");
      }

      
      /* Removes all the instances with weight equal to 0.
	 MUST be done since condition (8) of Keerthi's paper 
	 is made with the assertion Ci > 0 (See equation (3a). */
      Instances data = new Instances(insts, insts.numInstances());
      for(int i = 0; i < insts.numInstances(); i++){
	if(insts.instance(i).weight() > 0)
	  data.add(insts.instance(i));
      }
      if (data.numInstances() == 0) {
	throw new Exception("No training instances left after removing " + 
			    "instance with either a weight null or a missing class!");
      }
      insts = data;
      
    }

    m_onlyNumeric = true;
    if (!m_checksTurnedOff) {
      for (int i = 0; i < insts.numAttributes(); i++) {
	if (i != insts.classIndex()) {
	  if (!insts.attribute(i).isNumeric()) {
	    m_onlyNumeric = false;
	    break;
	  }
	}
      }
    }

    if (!m_checksTurnedOff) {
      m_Missing = new ReplaceMissingValues();
      m_Missing.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Missing); 
    } else {
      m_Missing = null;
    }

    if (!m_onlyNumeric) {
      m_NominalToBinary = new NominalToBinary();
      m_NominalToBinary.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_NominalToBinary);
    } else {
      m_NominalToBinary = null;
    }

    if (m_filterType == FILTER_STANDARDIZE) {
      m_Filter = new Standardize();
      m_Filter.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Filter); 
    } else if (m_filterType == FILTER_NORMALIZE) {
      m_Filter = new Normalize();
      m_Filter.setInputFormat(insts);
      insts = Filter.useFilter(insts, m_Filter); 
    } else {
      m_Filter = null;
    }

    m_classIndex = insts.classIndex();
    m_classAttribute = insts.classAttribute();

    // Generate subsets representing each class
    Instances[] subsets = new Instances[insts.numClasses()];
    for (int i = 0; i < insts.numClasses(); i++) {
      subsets[i] = new Instances(insts, insts.numInstances());
    }
    for (int j = 0; j < insts.numInstances(); j++) {
      Instance inst = insts.instance(j);
      subsets[(int)inst.classValue()].add(inst);
    }
    for (int i = 0; i < insts.numClasses(); i++) {
      subsets[i].compactify();
    }

    // Build the binary classifiers
    Random rand = new Random(m_randomSeed);
    m_classifiers = new BinarymySMO[insts.numClasses()][insts.numClasses()];
    for (int i = 0; i < insts.numClasses(); i++) {
      for (int j = i + 1; j < insts.numClasses(); j++) {
	m_classifiers[i][j] = new BinarymySMO();
	Instances data = new Instances(insts, insts.numInstances());
	for (int k = 0; k < subsets[i].numInstances(); k++) {
	  data.add(subsets[i].instance(k));
	}
	for (int k = 0; k < subsets[j].numInstances(); k++) {
	  data.add(subsets[j].instance(k));
	}
	data.compactify();
	data.randomize(rand);
	m_classifiers[i][j].buildClassifier(data, i, j, 
					    m_fitLogisticModels,
					    m_numFolds, m_randomSeed);
      }
    }
  }

  /**
   * Estimates class probabilities for given instance.
   */
  public double[] distributionForInstance(Instance inst) throws Exception {

    // Filter instance
    if (!m_checksTurnedOff) {
      m_Missing.input(inst);
      m_Missing.batchFinished();
      inst = m_Missing.output();
    }

    if (!m_onlyNumeric) {
      m_NominalToBinary.input(inst);
      m_NominalToBinary.batchFinished();
      inst = m_NominalToBinary.output();
    }
    
    if (m_Filter != null) {
      m_Filter.input(inst);
      m_Filter.batchFinished();
      inst = m_Filter.output();
    }
    
    if (!m_fitLogisticModels) {
      double[] result = new double[inst.numClasses()];
      for (int i = 0; i < inst.numClasses(); i++) {
	for (int j = i + 1; j < inst.numClasses(); j++) {
	  if ((m_classifiers[i][j].m_alpha != null) || 
	      (m_classifiers[i][j].m_sparseWeights != null)) {
	    double output = m_classifiers[i][j].mySVMOutput(-1, inst);
	    if (output > 0) {
	      result[j] += 1;
	    } else {
	      result[i] += 1;
	    }
	  }
	} 
      }
      Utils.normalize(result);
      return result;
    } else {

      // We only need to do pairwise coupling if there are more
      // then two classes.
      if (inst.numClasses() == 2) {
	double[] newInst = new double[2];
	newInst[0] = m_classifiers[0][1].mySVMOutput(-1, inst);
	newInst[1] = Instance.missingValue();
	return m_classifiers[0][1].m_logistic.
	  distributionForInstance(new Instance(1, newInst));
      }
      double[][] r = new double[inst.numClasses()][inst.numClasses()];
      double[][] n = new double[inst.numClasses()][inst.numClasses()];
      for (int i = 0; i < inst.numClasses(); i++) {
	for (int j = i + 1; j < inst.numClasses(); j++) {
	  if ((m_classifiers[i][j].m_alpha != null) || 
	      (m_classifiers[i][j].m_sparseWeights != null)) {
	    double[] newInst = new double[2];
	    newInst[0] = m_classifiers[i][j].mySVMOutput(-1, inst);
	    newInst[1] = Instance.missingValue();
	    r[i][j] = m_classifiers[i][j].m_logistic.
	      distributionForInstance(new Instance(1, newInst))[0];
	    n[i][j] = m_classifiers[i][j].m_sumOfWeights;
	  }
	}
      }
      return pairwiseCoupling(n, r);
    }
  }

  /**
   * Implements pairwise coupling.
   *
   * @param n the sum of weights used to train each model
   * @param r the probability estimate from each model
   * @return the coupled estimates
   */
  public double[] pairwiseCoupling(double[][] n, double[][] r) {

    // Initialize p and u array
    double[] p = new double[r.length];
    for (int i =0; i < p.length; i++) {
      p[i] = 1.0 / (double)p.length;
    }
    double[][] u = new double[r.length][r.length];
    for (int i = 0; i < r.length; i++) {
      for (int j = i + 1; j < r.length; j++) {
	u[i][j] = 0.5;
      }
    }

    // firstSum doesn't change
    double[] firstSum = new double[p.length];
    for (int i = 0; i < p.length; i++) {
      for (int j = i + 1; j < p.length; j++) {
	firstSum[i] += n[i][j] * r[i][j];
	firstSum[j] += n[i][j] * (1 - r[i][j]);
      }
    }

    // Iterate until convergence
    boolean changed;
    do {
      changed = false;
      double[] secondSum = new double[p.length];
      for (int i = 0; i < p.length; i++) {
	for (int j = i + 1; j < p.length; j++) {
	  secondSum[i] += n[i][j] * u[i][j];
	  secondSum[j] += n[i][j] * (1 - u[i][j]);
	}
      }
      for (int i = 0; i < p.length; i++) {
	if ((firstSum[i] == 0) || (secondSum[i] == 0)) {
	  if (p[i] > 0) {
	    changed = true;
	  }
	  p[i] = 0;
	} else {
	  double factor = firstSum[i] / secondSum[i];
	  double pOld = p[i];
	  p[i] *= factor;
	  if (Math.abs(pOld - p[i]) > 1.0e-3) {
	    changed = true;
	  }
	}
      }
      Utils.normalize(p);
      for (int i = 0; i < r.length; i++) {
	for (int j = i + 1; j < r.length; j++) {
	  u[i][j] = p[i] / (p[i] + p[j]);
	}
      }
    } while (changed);
    return p;
  }

  /**
   * Returns an array of votes for the given instance.
   * @param inst the instance
   * @return array of votex
   * @exception Exception if something goes wrong
   */
  public double obtainVotes(Instance inst) throws Exception {

    // Filter instance
    if (!m_checksTurnedOff) {
      m_Missing.input(inst);
      m_Missing.batchFinished();
      inst = m_Missing.output();
    }

    if (!m_onlyNumeric) {
      m_NominalToBinary.input(inst);
      m_NominalToBinary.batchFinished();
      inst = m_NominalToBinary.output();
    }
    
    if (m_Filter != null) {
      m_Filter.input(inst);
      m_Filter.batchFinished();
      inst = m_Filter.output();
    }

    double[] votes = new double[inst.numClasses()];
    for (int i = 0; i < inst.numClasses(); i++) {
      for (int j = i + 1; j < inst.numClasses(); j++) {
	double output = m_classifiers[i][j].mySVMOutput(-1, inst);
	if (output > 0) {
	  votes[j] += output;
	} else {
	  votes[i] += output;
	}
      }
    }
    double ss=votes[0] + votes[1];
    return ss;
  }

  /**
   * Returns the weights in sparse format.
   */
  public double [][][] sparseWeights() {
    
    int numValues = m_classAttribute.numValues();
    double [][][] sparseWeights = new double[numValues][numValues][];
    
    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	sparseWeights[i][j] = m_classifiers[i][j].m_sparseWeights;
      }
    }
    
    return sparseWeights;
  }
  
  /**
   * Returns the indices in sparse format.
   */
  public int [][][] sparseIndices() {
    
    int numValues = m_classAttribute.numValues();
    int [][][] sparseIndices = new int[numValues][numValues][];

    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	sparseIndices[i][j] = m_classifiers[i][j].m_sparseIndices;
      }
    }
    
    return sparseIndices;
  }
  
  /**
   * Returns the bias of each binary mySMO.
   */
  public double [][] bias() {
    
    int numValues = m_classAttribute.numValues();
    double [][] bias = new double[numValues][numValues];

    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	bias[i][j] = m_classifiers[i][j].m_b;
      }
    }
    
    return bias;
  }
  
  /*
   * Returns the number of values of the class attribute.
   */
  public int numClassAttributeValues() {

    return m_classAttribute.numValues();
  }
  
  /*
   * Returns the names of the class attributes.
   */
  public String [] classAttributeNames() {

    int numValues = m_classAttribute.numValues();
    
    String [] classAttributeNames = new String[numValues];
    
    for (int i = 0; i < numValues; i++) {
      classAttributeNames[i] = m_classAttribute.value(i);
    }
    
    return classAttributeNames;
  }
  
  /**
   * Returns the attribute names.
   */
  public String [][][] attributeNames() {
    
    int numValues = m_classAttribute.numValues();
    String [][][] attributeNames = new String[numValues][numValues][];
    
    for (int i = 0; i < numValues; i++) {
      for (int j = i + 1; j < numValues; j++) {
	int numAttributes = m_classifiers[i][j].m_data.numAttributes();
	String [] attrNames = new String[numAttributes];
	for (int k = 0; k < numAttributes; k++) {
	  attrNames[k] = m_classifiers[i][j].m_data.attribute(k).name();
	}
	attributeNames[i][j] = attrNames;          
      }
    }
    return attributeNames;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(13);

    newVector.addElement(new Option("\tThe complexity constant C. (default 1)",
				    "C", 1, "-C <double>"));
    newVector.addElement(new Option("\tThe exponent for the "
				    + "polynomial kernel. (default 1)",
				    "E", 1, "-E <double>"));
    newVector.addElement(new Option("\tGamma for the RBF kernel. (default 0.01)",
				    "G", 1, "-G <double>"));
    newVector.addElement(new Option("\tWhether to 0=normalize/1=standardize/2=neither. " +
				    "(default 0=normalize)",
				    "N", 1, "-N"));
    newVector.addElement(new Option("\tFeature-space normalization (only for\n" +
				    "\tnon-linear polynomial kernels).",
				    "F", 0, "-F"));
    newVector.addElement(new Option("\tUse lower-order terms (only for non-linear\n" +
				    "\tpolynomial kernels).",
				    "O", 0, "-O"));
    newVector.addElement(new Option("\tUse RBF kernel. " +
    				    "(default poly)",
				    "R", 0, "-R"));
    newVector.addElement(new Option("\tThe size of the kernel cache. " +
				    "(default 250007, use 0 for full cache)",
				    "A", 1, "-A <int>"));
    newVector.addElement(new Option("\tThe tolerance parameter. " +
				    "(default 1.0e-3)",
				    "L", 1, "-L <double>"));
    newVector.addElement(new Option("\tThe epsilon for round-off error. " +
				    "(default 1.0e-12)",
				    "P", 1, "-P <double>"));
    newVector.addElement(new Option("\tFit logistic models to SVM outputs. ",
				    "M", 0, "-M"));
    newVector.addElement(new Option("\tThe number of folds for the internal\n" +
				    "\tcross-validation. " +
				    "(default -1, use training data)",
				    "V", 1, "-V <double>"));
    newVector.addElement(new Option("\tThe random number seed. " +
				    "(default 1)",
				    "W", 1, "-W <double>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * -C num <br>
   * The complexity constant C. (default 1)<p>
   *
   * -E num <br>
   * The exponent for the polynomial kernel. (default 1) <p>
   *
   * -G num <br>
   * Gamma for the RBF kernel. (default 0.01) <p>
   *
   * -N <0|1|2> <br>
   * Whether to 0=normalize/1=standardize/2=neither. (default 0=normalize)<p>
   *
   * -F <br>
   * Feature-space normalization (only for  non-linear polynomial kernels). <p>
   *
   * -O <br>
   * Use lower-order terms (only for non-linear polynomial kernels). <p>
   *
   * -R <br>
   * Use RBF kernel (default poly). <p>
   * 
   * -A num <br> Sets the size of the kernel cache. Should be a prime
   * number. (default 250007, use 0 for full cache) <p>
   *
   * -L num <br>
   * Sets the tolerance parameter. (default 1.0e-3)<p>
   *
   * -P num <br> 
   * Sets the epsilon for round-off error. (default 1.0e-12)<p>
   *
   * -M <br>
   * Fit logistic models to SVM outputs.<p>
   *
   * -V num <br>
   * Number of folds for cross-validation used to generate data
   * for logistic models. (default -1, use training data)
   *
   * -W num <br>
   * Random number seed. (default 1)
   *
   * @param options the list of options as an array of strings
   * @exception Exception if an option is not supported 
   */
  public void setOptions(String[] options) throws Exception {
    
    String complexityString = Utils.getOption('C', options);
    if (complexityString.length() != 0) {
      m_C = (new Double(complexityString)).doubleValue();
    } else {
      m_C = 1.0;
    }
    String exponentsString = Utils.getOption('E', options);
    if (exponentsString.length() != 0) {
      m_exponent = (new Double(exponentsString)).doubleValue();
    } else {
      m_exponent = 1.0;
    }
    String gammaString = Utils.getOption('G', options);
    if (gammaString.length() != 0) {
      m_gamma = (new Double(gammaString)).doubleValue();
    } else {
      m_gamma = 0.01;
    }
    String cacheString = Utils.getOption('A', options);
    if (cacheString.length() != 0) {
      m_cacheSize = Integer.parseInt(cacheString);
    } else {
      m_cacheSize = 250007;
    }
    String toleranceString = Utils.getOption('L', options);
    if (toleranceString.length() != 0) {
      m_tol = (new Double(toleranceString)).doubleValue();
    } else {
      m_tol = 1.0e-3;
    }
    String epsilonString = Utils.getOption('P', options);
    if (epsilonString.length() != 0) {
      m_eps = (new Double(epsilonString)).doubleValue();
    } else {
      m_eps = 1.0e-12;
    }
    m_useRBF = Utils.getFlag('R', options);
    String nString = Utils.getOption('N', options);
    if (nString.length() != 0) {
      setFilterType(new SelectedTag(Integer.parseInt(nString), TAGS_FILTER));
    } else {
      setFilterType(new SelectedTag(FILTER_NORMALIZE, TAGS_FILTER));
    }
    m_featureSpaceNormalization = Utils.getFlag('F', options);
    if ((m_useRBF) && (m_featureSpaceNormalization)) {
      throw new Exception("RBF machine doesn't require feature-space normalization.");
    }
    if ((m_exponent == 1.0) && (m_featureSpaceNormalization)) {
      throw new Exception("Can't use feature-space normalization with linear machine.");
    }
    m_lowerOrder = Utils.getFlag('O', options);
    if ((m_useRBF) && (m_lowerOrder)) {
      throw new Exception("Can't use lower-order terms with RBF machine.");
    }
    if ((m_exponent == 1.0) && (m_lowerOrder)) {
      throw new Exception("Can't use lower-order terms with linear machine.");
    }
    m_fitLogisticModels = Utils.getFlag('M', options);
    String foldsString = Utils.getOption('V', options);
    if (foldsString.length() != 0) {
      m_numFolds = Integer.parseInt(foldsString);
    } else {
      m_numFolds = -1;
    }
    String randomSeedString = Utils.getOption('W', options);
    if (randomSeedString.length() != 0) {
      m_randomSeed = Integer.parseInt(randomSeedString);
    } else {
      m_randomSeed = 1;
    }
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [21];
    int current = 0;

    options[current++] = "-C"; options[current++] = "" + m_C;
    options[current++] = "-E"; options[current++] = "" + m_exponent;
    options[current++] = "-G"; options[current++] = "" + m_gamma;
    options[current++] = "-A"; options[current++] = "" + m_cacheSize;
    options[current++] = "-L"; options[current++] = "" + m_tol;
    options[current++] = "-P"; options[current++] = "" + m_eps;
    options[current++] = "-N"; options[current++] = "" + m_filterType;
    if (m_featureSpaceNormalization) {
      options[current++] = "-F";
    }
    if (m_lowerOrder) {
      options[current++] = "-O";
    }
    if (m_useRBF) {
      options[current++] = "-R";
    }
    if (m_fitLogisticModels) {
      options[current++] = "-M";
    }
    options[current++] = "-V"; options[current++] = "" + m_numFolds;
    options[current++] = "-W"; options[current++] = "" + m_randomSeed;    

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String exponentTipText() {
    return "The exponent for the polynomial kernel.";
  }
  
  /**
   * Get the value of exponent. 
   *
   * @return Value of exponent.
   */
  public double getExponent() {
    
    return m_exponent;
  }
  
  /**
   * Set the value of exponent. If linear kernel
   * is used, rescaling and lower-order terms are
   * turned off.
   *
   * @param v  Value to assign to exponent.
   */
  public void setExponent(double v) {
    
    if (v == 1.0) {
      m_featureSpaceNormalization = false;
      m_lowerOrder = false;
    }
    m_exponent = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String gammaTipText() {
    return "The value of the gamma parameter for RBF kernels.";
  }
  
  /**
   * Get the value of gamma. 
   *
   * @return Value of gamma.
   */
  public double getGamma() {
    
    return m_gamma;
  }
  
  /**
   * Set the value of gamma. 
   *
   * @param v  Value to assign to gamma.
   */
  public void setGamma(double v) {
    
    m_gamma = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String cTipText() {
    return "The complexity parameter C.";
  }
  
  /**
   * Get the value of C.
   *
   * @return Value of C.
   */
  public double getC() {
    
    return m_C;
  }
  
  /**
   * Set the value of C.
   *
   * @param v  Value to assign to C.
   */
  public void setC(double v) {
    
    m_C = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String toleranceParameterTipText() {
    return "The tolerance parameter (shouldn't be changed).";
  }
  
  /**
   * Get the value of tolerance parameter.
   * @return Value of tolerance parameter.
   */
  public double getToleranceParameter() {
    
    return m_tol;
  }
  
  /**
   * Set the value of tolerance parameter.
   * @param v  Value to assign to tolerance parameter.
   */
  public void setToleranceParameter(double v) {
    
    m_tol = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String epsilonTipText() {
    return "The epsilon for round-off error (shouldn't be changed).";
  }
  
  /**
   * Get the value of epsilon.
   * @return Value of epsilon.
   */
  public double getEpsilon() {
    
    return m_eps;
  }
  
  /**
   * Set the value of epsilon.
   * @param v  Value to assign to epsilon.
   */
  public void setEpsilon(double v) {
    
    m_eps = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String cacheSizeTipText() {
    return "The size of the kernel cache (should be a prime number). Use 0 for full cache.";
  }
  
  /**
   * Get the size of the kernel cache
   * @return Size of kernel cache.
   */
  public int getCacheSize() {
    
    return m_cacheSize;
  }
  
  /**
   * Set the value of the kernel cache.
   * @param v  Size of kernel cache.
   */
  public void setCacheSize(int v) {
    
    m_cacheSize = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String filterTypeTipText() {
    return "Determines how/if the data will be transformed.";
  }
  
  /**
   * Gets how the training data will be transformed. Will be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @return the filtering mode
   */
  public SelectedTag getFilterType() {

    return new SelectedTag(m_filterType, TAGS_FILTER);
  }
  
  /**
   * Sets how the training data will be transformed. Should be one of
   * FILTER_NORMALIZE, FILTER_STANDARDIZE, FILTER_NONE.
   *
   * @param newType the new filtering mode
   */
  public void setFilterType(SelectedTag newType) {
    
    if (newType.getTags() == TAGS_FILTER) {
      m_filterType = newType.getSelectedTag().getID();
    }
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useRBFTipText() {
    return "Whether to use an RBF kernel instead of a polynomial one.";
  }

  /**
   * Check if the RBF kernel is to be used.
   * @return true if RBF
   */
  public boolean getUseRBF() {
    
    return m_useRBF;
  }
  
  /**
   * Set if the RBF kernel is to be used.
   * @param v  true if RBF
   */
  public void setUseRBF(boolean v) {

    if (v) {
      m_featureSpaceNormalization = false;
      m_lowerOrder = false;
    }
    m_useRBF = v;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String featureSpaceNormalizationTipText() {
    return "Whether feature-space normalization is performed (only "
      + "available for non-linear polynomial kernels).";
  }
  
  /**
   * Check whether feature spaces is being normalized.
   * @return true if feature space is normalized.
   */
  public boolean getFeatureSpaceNormalization() throws Exception {

    return m_featureSpaceNormalization;
  }
  
  /**
   * Set whether feature space is normalized.
   * @param v  true if feature space is to be normalized.
   */
  public void setFeatureSpaceNormalization(boolean v) throws Exception {
    
    if ((m_useRBF) || (m_exponent == 1.0)) {
      m_featureSpaceNormalization = false;
    } else {
      m_featureSpaceNormalization = v;
    }
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String lowerOrderTermsTipText() {
    return "Whether lower order polyomials are also used (only "
      + "available for non-linear polynomial kernels).";
  }
  
  /**
   * Check whether lower-order terms are being used.
   * @return Value of lowerOrder.
   */
  public boolean getLowerOrderTerms() {
    
    return m_lowerOrder;
  }
  
  /**
   * Set whether lower-order terms are to be used. Defaults
   * to false if a linear machine is built.
   * @param v  Value to assign to lowerOrder.
   */
  public void setLowerOrderTerms(boolean v) {
    
    if (m_exponent == 1.0 || m_useRBF) {
      m_lowerOrder = false;
    } else {
      m_lowerOrder = v;
    }
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String buildLogisticModelsTipText() {
    return "Whether to fit logistic models to the outputs (for proper "
      + "probability estimates).";
  }

  /**
   * Get the value of buildLogisticModels.
   *
   * @return Value of buildLogisticModels.
   */
  public boolean getBuildLogisticModels() {
    
    return m_fitLogisticModels;
  }
  
  /**
   * Set the value of buildLogisticModels.
   *
   * @param newbuildLogisticModels Value to assign to buildLogisticModels.
   */
  public void setBuildLogisticModels(boolean newbuildLogisticModels) {
    
    m_fitLogisticModels = newbuildLogisticModels;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numFoldsTipText() {
    return "The number of folds for cross-validation used to generate "
      + "training data for logistic models (-1 means use training data).";
  }
  
  /**
   * Get the value of numFolds.
   *
   * @return Value of numFolds.
   */
  public int getNumFolds() {
    
    return m_numFolds;
  }
  
  /**
   * Set the value of numFolds.
   *
   * @param newnumFolds Value to assign to numFolds.
   */
  public void setNumFolds(int newnumFolds) {
    
    m_numFolds = newnumFolds;
  }
     
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {
    return "Random number seed for the cross-validation.";
  }
  
  /**
   * Get the value of randomSeed.
   *
   * @return Value of randomSeed.
   */
  public int getRandomSeed() {
    
    return m_randomSeed;
  }
  
  /**
   * Set the value of randomSeed.
   *
   * @param newrandomSeed Value to assign to randomSeed.
   */
  public void setRandomSeed(int newrandomSeed) {
    
    m_randomSeed = newrandomSeed;
  }
  
  /**
   * Prints out the classifier.
   *
   * @return a description of the classifier as a string
   */
  public String toString() {
    
    StringBuffer text = new StringBuffer();
    int printed = 0;
    
    if ((m_classAttribute == null)) {
      return "mySMO: No model built yet.";
    }
    try {
      text.append("mySMO\n\n");
      for (int i = 0; i < m_classAttribute.numValues(); i++) {
	for (int j = i + 1; j < m_classAttribute.numValues(); j++) {
	  text.append("Classifier for classes: " + 
		      m_classAttribute.value(i) + ", " +
		      m_classAttribute.value(j) + "\n\n");
	  text.append(m_classifiers[i][j]);
	  if (m_fitLogisticModels) {
	    text.append("\n\n");
	    if ( m_classifiers[i][j].m_logistic == null) {
	      text.append("No logistic model has been fit.\n");
	    } else {
	      text.append(m_classifiers[i][j].m_logistic);
	    }
	  }
	  text.append("\n\n");
	}
      }
    } catch (Exception e) {
      return "Can't print mySMO classifier.";
    }
    
    return text.toString();
  }

  public double myScore(double m_funcval,double m_sgnf)
  {
     double m_numt = 1.0;
     double m_expo = 1.0;
     double m_tfunc = 0.0;
     double m_fres = 0.0;
    
     if(m_sgnf > 0.0)
     {
        m_tfunc = m_tfunc - m_funcval;
     }else{
        m_tfunc = m_funcval;
     }

     if(m_funcval != 0.0)
     {
        m_expo = Math.exp(m_tfunc);
     }
     
     m_fres = m_numt / (m_numt + m_expo);

     return m_fres;
  }

  
  public void mygenstat(Instances m_tinst)
  {
     System.out.println("remove text till here...");
     String[] m_cnames = this.classAttributeNames();
     try{
      for (int m_k = 0; m_k < m_tinst.numInstances(); m_k++) {
	  Instance mc_inst = m_tinst.instance(m_k);
      String m_classname = mc_inst.stringValue(mc_inst.classAttribute());
     // Instance mc_wmiss = (Instance)mc_inst.copy();
     // mc_wmiss.setData(m_tinst);
     // double mpval = ((Classifier)classifier).classifyInstance(mc_wmiss);
      
      /*
      int m_rid=0;
            if(m_classname.equals(m_cnames[0]))
            {
               m_rid =0;
            }else if(m_classname.equals(m_cnames[1])){
               m_rid=1;
            }else{
               System.out.println("unknown class: m_classname " + m_classname + ", m_cname " + m_cnames[0] + ", " + m_cnames[1]);
               System.exit(1);
            }*/

//      System.out.println("classvalue " + mc_inst.classValue() + " , classname " + m_classname);
	  double m_output=obtainVotes(mc_inst);
	  double[][] m_bias=bias();
      double m_myscore = myScore(m_output,m_bias[0][1]);

     // if(mc_inst.classAttribute().isNumeric())
     // {
      //    System.out.println(m_k + "]class=" + m_classname + ",output=" + m_output + ",bias=" + m_bias[0][1] + ",score=" + m_myscore+",pred="+mpval);
      //}else{
      //    String m_pclass = mc_inst.classAttribute().value((int)mpval);
      //    double mprob=classifier.distributionForInstance(withMissing)[(int)mpval];
      //    System.out.println(m_k + "]class=" + m_classname + ",output=" + m_output + ",bias=" + m_bias[0][1] + ",score=" + m_myscore+",pred="+m_pclass+",prob="+m_prob);
      //}

      System.out.println(m_k + "]class=" + m_classname + ",output=" + m_output + ",bias=" + m_bias[0][1] + ",score=" + m_myscore);
      }
     }catch (Exception e) {
      e.printStackTrace();
      System.err.println(e.getMessage());
    }

  }
  
  public String[] m_myclasses;
  /**
   * Main method for testing this class.
   */
  public static void main(String[] argv) {

   //Classifier scheme;
   SMO scheme;

    try {
      scheme = new SMO();
     
         // String m_inpStr = "-C 1.0 -E 1.0 -G 0.01 -A 250007 -L 0.0010 -P 1.0E-12 -M -N 0 -V 10 -W 1 -t ../EXPERIMENT/filtered/2new/plasm-new-infogain50-dataset.arff -T ./plasm-overlap50-dataset.arff -p 0";
         String m_inpStr = "-C 1.0 -E 1.0 -G 0.01 -A 250007 -L 0.0010 -P 1.0E-12 -M -N 0 -V 10 -W 1 -i -t ./plasm-train184newinfo21-dataset.arff -T ./plasm-gametinfo21-dataset.arff";
//         String m_inpStr = "-C 1.0 -E 1.0 -G 0.01 -A 250007 -L 0.0010 -P 1.0E-12 -M -N 0 -R -V 10 -W 1 -t ../EXPERIMENT/filtered/2class/plasm-infogainranker-t60.arff -d ./smo.model";
//         String m_inpStr = "-T ./1test.arff -l ./smo.model";
         String[] m_res_inpStr = m_inpStr.split("\\s");
         System.out.println(Evaluation.evaluateModel(scheme,m_res_inpStr));

      
         m_res_inpStr = m_inpStr.split("\\s");
         String m_trainFileName = Utils.getOption('t',m_res_inpStr);
         String m_testFileName = Utils.getOption('T',m_res_inpStr);
         
         BufferedReader m_trainReader = null;
         BufferedReader m_testReader = null;
         
         Instances m_train = null;
         Instances m_test = null;

         if(m_trainFileName.length() != 0)
         {
            System.out.println("trainfile = " + m_trainFileName);
            m_trainReader = new BufferedReader(new FileReader(m_trainFileName));
            m_train = new Instances(m_trainReader);
            m_train.setClassIndex(m_train.numAttributes() - 1);
         }
         
          if(m_testFileName.length() != 0)
         {
            System.out.println("testfile = " + m_testFileName);
            m_testReader = new BufferedReader(new FileReader(m_testFileName));
            m_test = new Instances(m_testReader);
            m_test.setClassIndex(m_test.numAttributes() - 1);
         }
        

     // if(m_trainFileName.length() != 0)
      //  scheme.mygenstat(m_train);

      System.out.println("#################"); 
      
      if(m_testFileName.length() != 0)
        scheme.mygenstat(m_test);

      //System.out.println(Evaluation.evaluateModel(scheme, argv));

    } catch (Exception e) {
      e.printStackTrace();
      System.err.println(e.getMessage());
    }
  }

@Override
public double classifyInstance(Instance arg0) throws Exception {
	// TODO Auto-generated method stub
	return 0;
}
}