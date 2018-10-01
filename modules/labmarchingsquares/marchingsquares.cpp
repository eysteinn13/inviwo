/*********************************************************************
*  Author  : Himangshu Saikia, Wiebke Koepp, ...
*  Init    : Monday, September 11, 2017 - 12:58:42
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo
{

	// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
	const ProcessorInfo MarchingSquares::processorInfo_
	{
		"org.inviwo.MarchingSquares",      // Class identifier
		"Marching Squares",                // Display name
		"KTH Lab",                          // Category
		CodeState::Experimental,           // Code state
		Tags::None,                        // Tags
	};

	const ProcessorInfo MarchingSquares::getProcessorInfo() const
	{
		return processorInfo_;
	}

	struct Square
	{
		vec2 bottomLeft, topLeft, topRight, bottomRight;
	};

	double ** MarchingSquares::SampleFromGrid(int x, int y, int dim, const VolumeRAM * vr, const size3_t dims)
	{
		double ** res = new double *[dim];
		for (int x = 0; x < dim; x++)
			res[x] = new double[dim];
		
		int offset = dim / 2;
		int xFilterIdx = 0;		
		for (int i = x - offset; i <= x + offset; i++)
		{
			int yFilterIdx = 0;
			for (int j = y - offset; j <= y + offset; j++)
			{
				//Make sure our indices are witin our grid
				int xIdx = i;
				int yIdx = j;
				xIdx = xIdx < 0 ? 0 : xIdx;
				xIdx = xIdx > (dims.x - 1) ? (dims.x - 1) : xIdx;
				yIdx = yIdx < 0 ? 0 : yIdx;
				yIdx = yIdx >(dims.y - 1) ? (dims.y - 1) : yIdx;

				res[xFilterIdx][yFilterIdx] = getInputValue(vr, dims, xIdx, yIdx, false);

				yFilterIdx++;
			}
			xFilterIdx++;
		}

		return res;
	}

	// Function to create Gaussian filter 
	void MarchingSquares::FilterCreation(double GKernel[3][3])
	{
		// intialising standard deviation to 1.0 
		double sigma = propSigma;
		double r, s = 2.0 * sigma * sigma;

		// sum is for normalization 
		double sum = 0.0;

		// generating 5x5 kernel 
		for (int x = -1; x <= 1; x++) {
			for (int y = -1; y <= 1; y++) {
				r = sqrt(x * x + y * y);
				GKernel[x + 1][y + 1] = (exp(-(r * r) / s)) / (M_PI * s);
				sum += GKernel[x + 1][y + 1];
			}
		}

		// normalising the Kernel 
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				GKernel[i][j] /= sum;
	}

	MarchingSquares::MarchingSquares()
		:Processor()
		, inData("volumeIn")
		, meshOut("meshOut")
		, propShowGrid("showGrid", "Show Grid")
		, propDeciderType("deciderType", "Decider Type")
		, propMultiple("multiple", "Iso Levels")
		, propIsoValue("isovalue", "Iso Value")
		, propSigma("sigma", "Sigma")
		, propGaussFilter("gaussFilter", "Gauss Filter")
		, propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f),
			vec4(0.0f), vec4(1.0f), vec4(0.1f),
			InvalidationLevel::InvalidOutput, PropertySemantics::Color)
		, propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
		, propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
	{
		// Register ports
		addPort(inData);
		addPort(meshOut);

		// Register properties
		addProperty(propShowGrid);
		addProperty(propGridColor);
		addProperty(propGaussFilter);

		addProperty(propDeciderType);
		propDeciderType.addOption("midpoint", "Mid Point", 0);
		propDeciderType.addOption("asymptotic", "Asymptotic", 1);

		addProperty(propMultiple);

		propMultiple.addOption("single", "Single", 0);
		addProperty(propIsoValue);		
		addProperty(propIsoColor);

		addProperty(propSigma);
		propSigma.setMinValue(0.0);
		propSigma.setMaxValue(10);

		propMultiple.addOption("multiple", "Multiple", 1);
		addProperty(propNumContours);
		addProperty(propIsoTransferFunc);

		// The default transfer function has just two blue points
		propIsoTransferFunc.get().clearPoints();
		propIsoTransferFunc.get().addPoint(vec2(0.0f, 1.0f), vec4(1.0f, 0.0f, 0.0f, 1.0f));
		propIsoTransferFunc.get().addPoint(vec2(1.0f, 1.0f), vec4(0.0f, 0.0f, 1.0f, 1.0f));
		propIsoTransferFunc.setCurrentStateAsDefault();

		util::hide(propGridColor, propNumContours, propIsoTransferFunc);

		// Show the grid color property only if grid is actually displayed
		propShowGrid.onChange([this]()
		{
			if (propShowGrid.get())
			{
				util::show(propGridColor);
			}
			else
			{
				util::hide(propGridColor);
			}
		});

		// Show options based on display of one or multiple iso contours
		propMultiple.onChange([this]()
		{
			if (propMultiple.get() == 0)
			{
				util::show(propIsoValue, propIsoColor);
				util::hide(propNumContours, propIsoTransferFunc);
			}
			else
			{
				//util::hide(propIsoValue);
				//util::show(propIsoColor, propNumContours);

				// TODO (Bonus): Comment out above if you are using the transfer function
				// and comment in below instead
				util::hide(propIsoValue, propIsoColor);
				util::show(propNumContours, propIsoTransferFunc);
			}
		});

	}


	bool MarchingSquares::midpointDecider(size_t x, size_t y, const VolumeRAM * vr, const size3_t dims)
	{
		double bottomLeftVal = getInputValue(vr, dims, x, y, propGaussFilter.get());
		double topLeftVal = getInputValue(vr, dims, x, y + 1, propGaussFilter.get());
		double topRightVal = getInputValue(vr, dims, x + 1, y + 1, propGaussFilter.get());
		double bottomRightVal = getInputValue(vr, dims, x + 1, y, propGaussFilter.get());

		double centerVal = (bottomLeftVal + topLeftVal + topRightVal + bottomRightVal) / 4;

		return centerVal < propIsoValue;
	}

	double MarchingSquares::interpolate(float cord0, float cord1, float val0, float val1, double iso)
	{
		return ((iso * (cord1 - cord0)) - ((val0 * cord1) - (val1 * cord0))) / (val1 - val0);
	}


	void MarchingSquares::process()
	{
		if (!inData.hasData()) {
			return;
		}

		// This results in a shared pointer to a volume
		auto vol = inData.getData();

		// Extract the minimum and maximum value from the input data
		const double minValue = vol->dataMap_.valueRange[0];
		const double maxValue = vol->dataMap_.valueRange[1];

		// Set the range for the isovalue to that minimum and maximum
		propIsoValue.setMinValue(minValue);
		propIsoValue.setMaxValue(maxValue);

		// You can print to the Inviwo console with Log-commands:
		LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue << ".");
		// You can also inform about errors and warnings:
		// LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
		// LogProcessorError("I am letting you know about an error"); // Will print error message in red
		// (There is also LogNetwork...() and just Log...(), these display a different source,
		// LogProcessor...() for example displays the name of the processor in the workspace while
		// Log...() displays the identifier of the processor (thus with multiple processors of the
		// same kind you would not know which one the information is coming from

		// Retreive data in a form that we can access it
		const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
		const size3_t dims = vol->getDimensions();

		LogProcessorInfo("Number of vertices in the x direction " << dims.x);
		LogProcessorInfo("Number of vertices in the y direction " << dims.y);

		// Initialize mesh and vertices
		auto mesh = std::make_shared<BasicMesh>();
		std::vector<BasicMesh::Vertex> vertices;

		// Values within the input data are accessed by the function below
		// It's input is the VolumeRAM from above, the dimensions of the volume
		// and the indeces i and j of the position to be accessed where
		// i is in [0, dims.x-1] and j is in [0, dims.y-1]
		float valueat00 = getInputValue(vr, dims, 0, 0, false);
		//("Value at (0,0) is: " << valueat00);
		//LogProcessorInfo("Value at (1,3) is: " << getInputValue(vr, dims, 1, 3, false));
		//LogProcessorInfo("Value at (3,3) is: " << getInputValue(vr, dims, 3, 3, false));
		// You can assume that dims.z = 1 and do not need to consider others cases


		// TODO (Bonus) Gaussian filter
		// Our input is const, but you need to compute smoothed data and write it somewhere
		// Create an editable volume like this:
		// Volume volSmoothed(vol->getDimensions(), vol->getDataFormat());
		// auto vrSmoothed = volSmoothed.getEditableRepresentation<VolumeRAM>();
		// Values can be set with
		// vrSmoothed->setFromDouble(vec3(i,j,0), value);
		// getting values works with an editable volume as well
		// getInputValue(vrSmoothed, dims, 0, 0);

		if(propGaussFilter.get())
		{
			//Initialize our Gauss filter
			double gauss[3][3];
			FilterCreation(gauss);

			//Instantiate our result array 
			filtered = new double*[dims.x];
			for (int i = 0; i < dims.x; i++)
				filtered[i] = new double[dims.y];
			
			for (size_t x = 0; x < dims.x; x++)
			{
				for (size_t y = 0; y < dims.y; y++)
				{
					LogProcessorInfo("Value at " << x << ",  " << y << " is: " << getInputValue(vr, dims, x, y, false));
					double ** gridSample = SampleFromGrid(x, y, 3, vr, dims);
					double result = 0.0;
					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							double current = gridSample[i][j];
							result += gauss[i][j] * gridSample[i][j];
						}
					}
					for (int i = 0; i < 3; i++)
						delete[] gridSample[i];
					delete[] gridSample;

					filtered[x][y] = result;
				}
			}
		}
		// Grid

		// Properties are accessed with propertyName.get() 
		if (propShowGrid.get())
		{
			// TODO: Add grid lines of the given color 

			// The function drawLineSegments creates two vertices at the specified positions, 
			// that are placed into the Vertex vector defining our mesh. 
			// An index buffer specifies which of those vertices should be grouped into to make up lines/trianges/quads.
			// Here two vertices make up a line segment.

			auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
			for (size_t x = 0; x < dims.x; x++)
			{
				float coord = (float)x / (dims.x - 1);
				vec2 v1 = vec2(coord, 0.0);
				vec2 v2 = vec2(coord, 1.0);
				drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
			}
			for (size_t y = 0; y < dims.y; y++)
			{
				float coord = (float)y / (dims.y - 1);
				vec2 v1 = vec2(0.0, coord);
				vec2 v2 = vec2(1.0, coord);
				drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
			}
		}

		// Iso contours

		if (propMultiple.get() == 0)
		{
			// TODO: Draw a single isoline at the specified isovalue (propIsoValue) 
			// and color it with the specified color (propIsoColor)

			//Step 1: Find all pairs of adjacent vertices that have different signs
			for (size_t xRaw = 0; xRaw < dims.x - 1; xRaw++)
			{
				for (size_t yRaw = 0; yRaw < dims.y - 1; yRaw++)
				{
					float x0 = (float)xRaw / (float)(dims.x - 1);
					float x1 = ((float)xRaw + 1) / (float)(dims.x - 1);
					float y0 = (float)yRaw / (float)(dims.y - 1);
					float y1 = ((float)yRaw + 1) / (float)(dims.y - 1);
					auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
					//Consider we are looking at squares
					vec3 bottomLeft(x0, y0, getInputValue(vr, dims, xRaw, yRaw, propGaussFilter.get()) < propIsoValue ? 0 : 1);
					vec3 topLeft(x0, y1, getInputValue(vr, dims, xRaw, yRaw + 1, propGaussFilter.get()) < propIsoValue ? 0 : 1);
					vec3 topRight(x1, y1, getInputValue(vr, dims, xRaw + 1, yRaw + 1, propGaussFilter.get()) < propIsoValue ? 0 : 1);
					vec3 bottomRight(x1, y0, getInputValue(vr, dims, xRaw + 1, yRaw, propGaussFilter.get()) < propIsoValue ? 0 : 1);
					double bottomRightVal = getInputValue(vr, dims, xRaw + 1, yRaw, propGaussFilter.get());
					double topRightVal = getInputValue(vr, dims, xRaw + 1, yRaw + 1, propGaussFilter.get());
					double bottomLeftVal = getInputValue(vr, dims, xRaw, yRaw, propGaussFilter.get());
					double topLeftVal = getInputValue(vr, dims, xRaw, yRaw + 1, propGaussFilter.get());
					double blbr = interpolate(x0, x1, bottomLeftVal, bottomRightVal, propIsoValue);
					double bltl = interpolate(y0, y1, bottomLeftVal, topLeftVal, propIsoValue);
					double trbr = interpolate(y0, y1, bottomRightVal, topRightVal, propIsoValue);
					double trtl = interpolate(x0, x1, topLeftVal, topRightVal, propIsoValue);


					LogProcessorInfo("Z vals " << topLeft.z << ", " << topRight.z << ", " << bottomRight.z << ", " << bottomLeft.z);

					//No nodes differ
					if (bottomLeft.z == topLeft.z && topLeft.z == topRight.z && topRight.z == bottomRight.z) {
						//Do nothing
					}

					//One node differs
					else if (bottomLeft.z == topLeft.z && topLeft.z == topRight.z) {
						drawLineSegment(vec2(blbr, y0), vec2(x1, trbr), propIsoColor.get(), indexBufferGrid, vertices);
					}
					else if (bottomLeft.z == topLeft.z && topLeft.z == bottomRight.z) {
						drawLineSegment(vec2(trtl, y1), vec2(x1, trbr), propIsoColor.get(), indexBufferGrid, vertices);
					}

					else if (bottomLeft.z == topRight.z && topRight.z == bottomRight.z) {
						drawLineSegment(vec2(trtl, y1), vec2(x0, bltl), propIsoColor.get(), indexBufferGrid, vertices);

					}
					else if (topLeft.z == topRight.z && topRight.z == bottomRight.z) {
						drawLineSegment(vec2(blbr, y0), vec2(x0, bltl), propIsoColor.get(), indexBufferGrid, vertices);
					}

					//Two nodes differ, case 1
					else if (topLeft.z == bottomLeft.z) {
						drawLineSegment(vec2(trtl, y1), vec2(blbr, y0), propIsoColor.get(), indexBufferGrid, vertices);
					}
					else if (topLeft.z == topRight.z) {
						drawLineSegment(vec2(x1, trbr), vec2(x0, bltl), propIsoColor.get(), indexBufferGrid, vertices);
					}

					//Two nodes differ, case 2== 0)
					else if (bottomLeft.z == topRight.z) {
						bool underIso = midpointDecider(xRaw, yRaw, vr, dims);
						if (propDeciderType.get() == 0)
						{
							if ((underIso && bottomLeft.z == 0) || (!underIso && topLeft.z == 0))
							{
								//Cut topLeft and BottomRight
								drawLineSegment(vec2(blbr, y0), vec2(x1, trbr), propIsoColor.get(), indexBufferGrid, vertices);
								drawLineSegment(vec2(trtl, y1), vec2(x0, bltl), propIsoColor.get(), indexBufferGrid, vertices);
							}
							else
							{
								//Cut bottomLeft and topRight
								drawLineSegment(vec2(blbr, y0), vec2(x0, bltl), propIsoColor.get(), indexBufferGrid, vertices);
								drawLineSegment(vec2(trtl, y1), vec2(x1, trbr), propIsoColor.get(), indexBufferGrid, vertices);

							}
						}
						else if (propDeciderType.get() == 1)
						{
							double points[4] = { blbr, trtl, x0, x1 };
							std::sort(points, points + 4);
							double correctSecond = points[1] == blbr ? y0 : y1;
							double correctThird = points[2] == trtl ? y1 : y0;
							drawLineSegment(vec2(points[0], bltl), vec2(points[1], correctSecond), propIsoColor.get(), indexBufferGrid, vertices);
							drawLineSegment(vec2(points[2], correctThird), vec2(points[3], trbr), propIsoColor.get(), indexBufferGrid, vertices);
						}
					}
				}
			}

		}
		else
		{
			// TODO: Draw the given number (propNumContours) of isolines between 
			// the minimum and maximum value
			double range = maxValue - minValue;
			double intervalSize = range / (propNumContours + 1);
			double contourValue = minValue + intervalSize;

			LogProcessorInfo("Range of values: " << range);
			LogProcessorInfo("Distance between contour lines: " << intervalSize);

			while (contourValue < maxValue)
			{
				LogProcessorInfo("Contour value: " << contourValue);
				/*--------------------------------------*/

				vec4 color = propIsoTransferFunc.get().sample(contourValue / maxValue);
				for (size_t xRaw = 0; xRaw < dims.x - 1; xRaw++)
				{
					for (size_t yRaw = 0; yRaw < dims.y - 1; yRaw++)
					{
						float x0 = (float)xRaw / (float)(dims.x - 1);
						float x1 = ((float)xRaw + 1) / (float)(dims.x - 1);
						float y0 = (float)yRaw / (float)(dims.y - 1);
						float y1 = ((float)yRaw + 1) / (float)(dims.y - 1);
						auto indexBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

						//Consider we are looking at squares
						vec3 bottomLeft(x0, y0, getInputValue(vr, dims, xRaw, yRaw, propGaussFilter.get()) < contourValue ? 0 : 1);
						vec3 topLeft(x0, y1, getInputValue(vr, dims, xRaw, yRaw + 1, propGaussFilter.get()) < contourValue ? 0 : 1);
						vec3 topRight(x1, y1, getInputValue(vr, dims, xRaw + 1, yRaw + 1, propGaussFilter.get()) < contourValue ? 0 : 1);
						vec3 bottomRight(x1, y0, getInputValue(vr, dims, xRaw + 1, yRaw, propGaussFilter.get()) < contourValue ? 0 : 1);

						double bottomRightVal = getInputValue(vr, dims, xRaw + 1, yRaw, propGaussFilter.get());
						double topRightVal = getInputValue(vr, dims, xRaw + 1, yRaw + 1, propGaussFilter.get());
						double bottomLeftVal = getInputValue(vr, dims, xRaw, yRaw, propGaussFilter.get());
						double topLeftVal = getInputValue(vr, dims, xRaw, yRaw + 1, propGaussFilter.get());

						double blbr = interpolate(x0, x1, bottomLeftVal, bottomRightVal, contourValue);
						double bltl = interpolate(y0, y1, bottomLeftVal, topLeftVal, contourValue);
						double trbr = interpolate(y0, y1, bottomRightVal, topRightVal, contourValue);
						double trtl = interpolate(x0, x1, topLeftVal, topRightVal, contourValue);

						//LogProcessorInfo("Z vals " << topLeft.z << ", " << topRight.z << ", " << bottomRight.z << ", " << bottomLeft.z);

						//No nodes differ
						if (bottomLeft.z == topLeft.z && topLeft.z == topRight.z && topRight.z == bottomRight.z) {
							//Do nothing
						}

						//One node differs
						else if (bottomLeft.z == topLeft.z && topLeft.z == topRight.z) {
							drawLineSegment(vec2(blbr, y0), vec2(x1, trbr), color, indexBufferGrid, vertices);
						}
						else if (bottomLeft.z == topLeft.z && topLeft.z == bottomRight.z) {
							drawLineSegment(vec2(trtl, y1), vec2(x1, trbr), color, indexBufferGrid, vertices);
						}

						else if (bottomLeft.z == topRight.z && topRight.z == bottomRight.z) {
							drawLineSegment(vec2(trtl, y1), vec2(x0, bltl), color, indexBufferGrid, vertices);

						}
						else if (topLeft.z == topRight.z && topRight.z == bottomRight.z) {
							drawLineSegment(vec2(blbr, y0), vec2(x0, bltl), color, indexBufferGrid, vertices);
						}

						//Two nodes differ, case 1
						else if (topLeft.z == bottomLeft.z) {
							drawLineSegment(vec2(trtl, y1), vec2(blbr, y0), color, indexBufferGrid, vertices);
						}
						else if (topLeft.z == topRight.z) {
							drawLineSegment(vec2(x1, trbr), vec2(x0, bltl), color, indexBufferGrid, vertices);
						}

						//Two nodes differ, case 2== 0)
						else if (bottomLeft.z == topRight.z) {
							bool underIso = midpointDecider(xRaw, yRaw, vr, dims);
							if (propDeciderType.get() == 0)
							{
								if ((underIso && bottomLeft.z == 0) || (!underIso && topLeft.z == 0))
								{
									//Cut topLeft and BottomRight
									drawLineSegment(vec2(blbr, y0), vec2(x1, trbr), color, indexBufferGrid, vertices);
									drawLineSegment(vec2(trtl, y1), vec2(x0, bltl), color, indexBufferGrid, vertices);
								}
								else
								{
									//Cut bottomLeft and topRight
									drawLineSegment(vec2(blbr, y0), vec2(x0, bltl), color, indexBufferGrid, vertices);
									drawLineSegment(vec2(trtl, y1), vec2(x1, trbr), color, indexBufferGrid, vertices);

								}
							}
							else if (propDeciderType.get() == 1)
							{
								double points[4] = { blbr, trtl, x0, x1 };
								std::sort(points, points + 4);
								double correctSecond = points[1] == blbr ? y0 : y1;
								double correctThird = points[2] == trtl ? y1 : y0;
								drawLineSegment(vec2(points[0], bltl), vec2(points[1], correctSecond), color, indexBufferGrid, vertices);
								drawLineSegment(vec2(points[2], correctThird), vec2(points[3], trbr), color, indexBufferGrid, vertices);
							}
						}
					}
				}
				contourValue += intervalSize;

				/*--------------------------------------*/
			}

			// TODO (Bonus): Use the transfer function property to assign a color
			// The transfer function normalizes the input data and sampling colors
			// from the transfer function assumes normalized input, that means
			// vec4 color = propIsoTransferFunc.get().sample(0.0f);
			// is the color for the minimum value in the data
			// vec4 color = propIsoTransferFunc.get().sample(1.0f);
			// is the color for the maximum value in the data


		}

		// Note: It is possible to add multiple index buffers to the same mesh,
		// thus you could for example add one for the grid lines and one for
		// each isoline
		// Also, consider to write helper functions to avoid code duplication
		// e.g. for the computation of a single iso contour

		mesh->addVertices(vertices);
		meshOut.setData(mesh);
	}

	double MarchingSquares::getInputValue(const VolumeRAM* data, const size3_t dims,
		const size_t i, const size_t j, bool filter = false) {
		
		if (filter)
		{
			return filtered[i][j];
		}

		// Check if the indices are withing the dimensions of the volume
		if (i < dims.x && j < dims.y) {
			return data->getAsDouble(size3_t(i, j, 0));
		}
		else {
			LogProcessorError(
				"Attempting to access data outside the boundaries of the volume, value is set to 0");
			LogProcessorError("Values tried: " << i << j);
			return 0;
		}
	}

	void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
		IndexBufferRAM* indexBuffer,
		std::vector<BasicMesh::Vertex>& vertices) {
		// Add first vertex
		indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
		// A vertex has a position, a normal, a texture coordinate and a color
		// we do not use normal or texture coordinate, but still have to specify them
		vertices.push_back({ vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color });
		// Add second vertex
		indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
		vertices.push_back({ vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color });
	}

} // namespace
