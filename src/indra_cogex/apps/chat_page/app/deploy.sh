#!/bin/bash

S3_URI_BASE=s3://bigmech
PATH_DEFAULT="/chat"

helpFunction() {

  echo ""
  echo "Usage: $0 [-i INDRALAB_VUE_TGZ] [-s S3_PATH]"
  echo -e "\t-i Path to indralab-vue tgz file. If provided, it will be installed. If not, it is assumed to already be installed."
  echo -e "\t-s The path in the bigmech S3 bucket to deploy to. Default is ${PATH_DEFAULT}. Path must start with /."
  echo -e "\t-h Help"
  echo ""
  echo "Example: $0 -i indralab-vue-0.1.0.tgz -s /chat"
  echo ""
  echo "This script will deploy the chat page to the bigmech S3 bucket. It assumes that a local packaging (*.tgz file)"
  echo "of indralab-vue is available for install or is already installed. This script assumes Node.js and the AWS CLI"
  echo "are installed and available."
  exit 1 # Exit script after printing help
}

while getopts "i:s:h" opt; do
  case "$opt" in
  i) indralabVueTgz="$OPTARG" ;;
  s) s3Path="$OPTARG" ;;
  h) helpFunction ;;
  ?) helpFunction ;; # Print helpFunction in case parameter is non-existent
  esac
done

# Print helpFunction in case parameters are empty (not needed for optional parameters)
#if [ -z "$indralabVueTgz" ] || [ -z "$s3Path" ]
#then
#   echo "Some or all of the parameters are empty";
#   helpFunction
#fi

# Install indralab-vue if provided
if [ -n "$indralabVueTgz" ]; then
  echo "Installing indralab-vue from ${indralabVueTgz}"
  npm install "${indralabVueTgz}"
fi

# Run the regular install
echo "Installing npm dependencies"
npm install

# Build the app
echo "Building index.html for chat app and dist"
npm run build

# Deploy to S3
if [ -n "$s3Path" ]; then
  PATH="$s3Path"
else # Use the default S3 URI
  PATH=$PATH_DEFAULT
fi

S3_URI="${S3_URI_BASE}${PATH}"
echo "Deploying to $S3_URI"

# Copy the content of the dist directory to the S3 bucket
# See https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3/sync.html for more details
aws s3 sync --exact-timestamps --delete dist/ "${S3_URI}"
echo "Deployment complete"
echo ""

# Force a cache refresh for CloudFront. NOTE: This requires the SUDO role
#aws cloudfront create-invalidation --distribution-id EFROMZ1D89URP --paths "${PATH}*"

# Instructions for manual cache invalidation
echo "Invalidate the CloudFront cache manually by going to https://us-east-1.console.aws.amazon.com/cloudfront/v3/home?region=us-east-1#/distributions and then:"
echo "  - Log in to the AWS console if needed"
echo "  - Click the distribution for discovery.indra.bio"
echo "  - Click the invalidations tab"
echo "  - Click 'Create invalidation'"
echo "  - Add the path /chat* to the list of paths to be invalidated"
echo "  - Click 'Create invalidation'"
echo "  - Wait for the invalidation to complete"
echo ""
