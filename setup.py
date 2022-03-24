import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wavelets_CMI-andreu-andreu",
    version="0.0.2",
    author="andreu andreu",
    author_email="author@example.com",
    description="wavelet decomposition and transfer entropy tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/andreuandreu/wavelets_CMI",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)