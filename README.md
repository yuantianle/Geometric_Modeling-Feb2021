# Geometric_Modeling-Feb2021

## [Project 1 - Mesh Subdivision & 3D Check Board Texturing](https://github.com/yuantianle/Geometric_Modeling-Feb2021/tree/main/Project1_MeshSubdivision_%26_3DCheckBoardTexturing)

### 1.1 Mesh Subdivision

|10k Tris|40k Tris|160k Tris|
|-|-|-|
|<img src="https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e8907b79-2980-41d1-9bf2-e89fed401474" width="420"/>|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/54ba6186-479f-45ab-b50b-78ffc8b2a7a4)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/9b1c29b7-1d28-464a-9a75-1c3f8aa211db)|

### 1.2 3D Check Board Texturing

|L = 1 | L = 0.5 | L = 0.25 |
|-|-|-|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/35e83bca-161a-4a82-a79c-7ac1ec46781f)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/601e2198-7f2d-4509-ba29-25409d5d223e)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/9f49be42-7033-47f7-9d5c-9d132ed30010)|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/7859b2b8-a297-4190-a2ab-8f73a2839c71)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e1bc791f-1067-445e-bd6b-f734c6abe27c)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d2cdc0be-cea0-42b0-a2f5-3b5d2e6282a1)|

## [Project 2 - Mesh Smoothing & 3D Heat Diffusion & 3D Texture Synthesis](https://github.com/yuantianle/Geometric_Modeling-Feb2021/tree/main/Project2_MeshSmoothing_%26_3DHeatDiffusion_%26_3DTexSynthesis)

### 2.1 Mesh Smoothing

<img src="https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/2a918247-4abc-4eb2-8193-7da494e1a985" width="180"/>

|Method|Angle View| Side View| Front View| Bottom View|
|-|-|-|-|-|
|Uniform scheme|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/92c2b9ea-bb9b-4f51-9fe0-813deb566b71)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d6ca4654-2847-42ea-9595-8f169428438b)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/6ddbec13-5ce8-4391-9e87-24d9c74a4f80)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d1a6fe8e-55ed-4ff0-ad72-d86e452544c9)|
|Cord scheme|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/f3ba182f-22d0-4a38-946d-da9293728448)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0232c0d7-627f-43f6-9d2c-0b360f6373cc)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/566c8d9e-f300-4a1a-a9f9-9dc47ef377fd)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/8dd2997b-5b65-4ba1-bd35-569baf63dc18)|
|Mean Curvature Flow scheme|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/42c76a98-478d-416a-97f1-ef43d3d8c038)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/bb386318-c2b6-4cef-a25e-dfaadeabaadb)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/375f26a2-e530-4648-881b-6d9824e1e16d)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/81dc3b25-28ed-4b16-8631-8fe4e6692406)|
|Mean Value Coordinates scheme|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c96691d4-4b6b-47d4-89ae-513f20b1b32f)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/f20a977f-80f2-4401-baaf-b7df89c222b8)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e4095b68-725b-4ad7-8a90-71f3b2f56838)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c7b73c99-17de-46da-b4c3-d0308dde8452)|

### 2.2 3D Heat Diffusion

|Iterating Solver|Result|Time cost|
|-|-|-|
|Gauss-Seidel (10000 times)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0875be6a-9ade-4a33-8883-953351b428bf)|457.678s|
|Explicit (10000 times)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/ec424c7b-f9f4-4f9c-a5da-ba4524885f10)|442.608s|
|Implicit (10000 times)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d51c939f-b9af-4643-8af1-f84a82fab109)|7.67s|

### 2.3 3D Texture Synthesis

|Texture|Left-Side View|Front View|Right-Side View|
|-|-|-|-|
|<img src="https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/7910d5f7-a8be-4684-9322-5b8d177dfbb7" width="45"/>|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/1b735869-03d1-4edb-ae75-8cf8bf382228)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/4f48b958-6756-4833-8e4b-d6580f7cf828)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/179afbb2-1fa4-471b-ab45-7e438c3f8f9d)|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0fa2c691-6247-4e10-8b57-1569a92b0e60)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d812307a-498c-406e-afbf-972c2777ba3b)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/2308ae6a-0b50-4f13-944b-4abd8f020c20)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/502e5539-e15f-40d8-a098-68b0cdc395e3)
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/eb7ab807-fcce-479c-8bf6-1530b8d40541)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/66f8bf80-9334-432a-b1b0-cddb1cce6a02)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/a7424f88-ec0b-46ca-bd76-26c9e1848c15)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/ef118007-dade-4f1b-8cbb-a61f7fb9dea6)|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/60eed796-983d-4914-b52a-6cdb14fa5a5d)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/b2e47195-82e7-4a1a-ac22-686796496d29)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/ae9c32ae-2a22-4983-a5e7-feefe972e396)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/649d34a8-adbf-4280-a73e-57e492671c1d)|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/07a7e501-dd7c-497c-864c-918c8338ad4f)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0d061203-4bc2-41de-bea4-5094923b0d15)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e19afa3e-a161-4efd-aeb1-9adabc5b2686)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/b4cef969-7902-46ee-bd83-0fbfa4c09593)|
|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/6db6c7ab-f8fe-463b-9b93-c747f369d95f)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c77e798c-b5d2-4387-ab9a-3c3f6fe07657)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/948aed73-fd52-4220-afa5-2582c9239399)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/9d03a61e-365a-402e-be22-7ba66b6cb641)|


## [Project 3 - Sketch Stylization & Curvature Tensor](https://github.com/yuantianle/Geometric_Modeling-Feb2021/tree/main/Project3_SketchStylization_%26_CurvatureTensor)

### 3.1 Curvature Tensor Visualization 

|Curvature Type|Donut Surface|Dragon Surface|
|-|-|-|
|Mean curvature|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/69e73123-e686-44a0-9e1c-be6003f2e899)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/f481709a-3a0e-4dec-8c3f-1568bea40ef0)|
|Gaussian curvature|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e8cdba9f-a414-4c61-8e97-80858c04a9f8)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0ea65afb-6350-4c6c-8425-e517e007be5c)|

#### Before Smoothing 

|Eigen Vector Type|Donut Surface|Dragon Surface|
|-|-|-|
|Minor principle|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/6b148f09-4f8f-43bf-8303-0fce052846b0)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/30a44a60-dbdc-4806-8662-b4830cd136dd)|
|Major principle|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/56a4440f-c0b5-4d76-85ae-3bb7381af983)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/09b22beb-4fb4-4a64-9a5d-20b144bb786b)|
|Duplicate principles|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/4bba3553-9480-4254-9804-a6e193fad01f)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e8d78b57-e005-4e55-b3ca-cb65b10aa398)|

#### After Smoothing 

|Eigen Vector Type|Donut Surface|Dragon Surface|
|-|-|-|
|Minor principle|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/d7c0d71d-5bca-45c2-baa8-f47b7c3a1e6d)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/3eb8bad5-0487-4d5e-ae1c-9b017b88a8ce)|
|Major principle|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/819f2cd3-acdd-4a09-910e-0f9b4e853d43)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/dbc90056-3c58-4645-8bc5-b8853d210349)|
|Duplicate principles|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c8da5fc9-a8a5-4f9c-b793-4a8bc2aaede8)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/2c7f85bb-cef3-4234-8d3d-da640b73f6f6)|

### 3.1 Sketch Stylization

|Mesh Type|Curvature Regression Line|Add Silhouette Stroke| Turn on Light |
|-|-|-|-|
|Donut|<img src="https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/ea52931f-8231-4e2f-b427-76a7024dc3aa" width="220"/>|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/1c27c6a0-1580-4cbc-80d2-c7f468e4bcf4)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/84b935de-57ed-4b64-9536-9915e7788748)|
|Bunny|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/6d59a3c1-28f0-43fa-8d1b-7bd1909dc6e0)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e236b198-70e8-4b3b-a78c-4b1203c51407)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/0bf8634f-6a80-4cdb-81bd-9d0a2cc35844)|
|Feline|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/1fecc8b0-a137-4a12-a819-de3104e11925)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/41783007-61cc-4fda-967b-1bb9362e6058)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/965b4ef6-81d7-40e7-82ac-4046ee003539)|
|Dragon|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/246e39bd-6d27-49a4-86b8-24989dd30500)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/9dd75048-303a-47a5-a999-962d7c75838a)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/560ccf0c-20c2-4076-ba46-a02d12c6d693)|
|Budda|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/691962ce-82dc-4452-81c4-171f701b36ca)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/7333a275-9bd3-4bba-8f0a-a0cffe7d2117)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c4565996-7b8c-4f75-bd5f-5fbe02b54c3f)|

## [Project 4 - Linear Rotation Invariant for Meshes Deformation](https://github.com/yuantianle/Geometric_Modeling-Feb2021/tree/main/Project4_LinearRotationInvariant_%26_MeshesDeformation)

### Transformation Results
|Transformation|Bunny Before|Bunny After|
|-|-|-|
|Rotate|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/1f1aee11-0f0b-4727-9d67-2bb8b8b6d1d0)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/f7d6a94a-8a8e-4b80-bd6c-7cd2fcdedebb)|
|Scale|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/eb42002a-df7e-45b1-9424-5e2344ab402c)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/e30e5ae2-0e23-4708-9921-17f8e60351e7)|
|Translate|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/2dcb1d87-256a-47f4-b7a4-ae90a96f9cfe)|![image](https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/c117a6c7-8e7d-4672-82ca-37c94c912763)|

<img src="https://github.com/yuantianle/Geometric_Modeling-Feb2021/assets/61530469/df838260-1eb6-41b5-978d-98c5031aa901"  width="90%"/>
