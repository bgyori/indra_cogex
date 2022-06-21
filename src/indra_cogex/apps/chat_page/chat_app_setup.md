### AWS S3, Cloudfront setup

A somewhat outdated, but useful video [here](https://www.youtube.com/watch?v=5r9Q-tI7mMw).
Helpful [SO question](https://stackoverflow.com/questions/31017105/how-do-you-set-a-default-root-object-for-subdirectories-for-a-statically-hosted), 
especially [this](https://stackoverflow.com/a/69157535/10478812) answer.

#### Build dist for indralab-vue and then for the chat app 

The full build process is documented in the Vue [readme file](./app/README.md).

To summarize:
- Build the `indralab-vue` dist by running `npm run build` in the root directory of the `indralab-vue` repository
- Pack the build by running `npm pack`. This will create a .tgz file that can be used as a standalone local installation.
- Install the `indralab-vue` dist: in the root directory of the Vue chat app, run `npm install /path/to/indralab-vue/indralab-vue-0.0.1.tgz`
- Build the Vue chat app: in the root directory of the Vue chat app, run `npm run build`

#### Upload app to S3

- Upload the _content_ of the `dist` directory for the chat app into a directory on S3 **with the same name as is 
  configured as `publicPath` in the Vue app's `vue.config.js` file**.
- Make sure that the content is publicly accessible (and that the bucket allows public access).

#### Set up a new origin in CloudFront for the content served from the S3 bucket

In the Distributions console for `discovery.indra.bio`, go to the `Origins` tab and click "Create origin":

- Click on "Origin domain". A list with AWS resources will appear. Find and select the S3 bucket with the uploaded content.
- _Leave the "Origin path" field empty!_
- Give the origin a name.
- Click "Create origin" and wait for the origin to be created.

#### Associate the new origin with a new Behavior

In the Distribution console, click on the **Behaviors** tab, then click on the **Create behavior** button.

- For path-patten, use the same name as the subfolder of the S3 bucket with a `*` added in the end, here that would 
  be `discovery_chat_app*`. This will allow matches for `discovery_chat_app`, `discovery_chat_app/` 
  and`discovery_chat_app/index.html`.
- For "origin and origin groups", select the S3 bucket origin created above, here that would be `cogex chat app`.
- For allowed HTTP methods, select `GET, HEAD`, since the chat app is only served via GET.
- Since the bucket allows public access, select `No` for "Restrict viewer access". _To restrict viewer access for a 
  private bucket, more setup is needed; see more
  [here](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/private-content-restricting-access-to-s3.html)._

#### Create a CloudFront Function to rewrite the URI

See the AWS documentation
[here](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/functions-tutorial.html) and
[here](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/functions-tutorial.html).

By using a CloudFront Function it's possible to rewrite requests that go to `/subfolder/` or `/subfolder` to be 
redirected to the URI `/subfolder/index.html` that's needed to load the content on S3. This makes it a lot easier to 
generate links to the content.

Steps:
- Go to the CloudFront console
- Click "Functions" in the left-hand menu (three stripes)
- If there isn't a function already, click "Create Function"
- Enter name and description
- Click "Create Function"
- For the JavaScript code, copy the code below and paste it into the editor

```javascript
function handler(event) {
    // NOTE: This example function is for a viewer request event trigger.
    // Choose viewer request for event trigger when you associate this function with a distribution.
    var request = event.request;
    var uri = request.uri;
    
    // Check whether the URI is missing a file name.
    if (uri.endsWith('/')) {
        request.uri += 'index.html';
    } 
    // Check whether the URI is missing a file extension.
    else if (!uri.includes('.')) {
        request.uri += '/index.html';
    }

    return request;
}
```

The code above is inspired by 
[this](https://github.com/aws-samples/amazon-cloudfront-functions/tree/main/url-rewrite-single-page-apps) repository.
To _redirect_ with 301 status from e.g. `/subfolder` to `/subfolder/` (so that it shows up in the browser) instead of 
_rewriting_ the request URI, see this
[link](https://github.com/digital-sailors/standard-redirects-for-cloudfront). Note that the code there is Node.js, 
not JavaScript and also is intended to run on Lambda, not CloudFront Functions, but the logic that's implemented there 
would be the same.

For full documentation on the event structure, see the
[AWS documentation](https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/functions-event-structure.html).

Once done with entering the code:
- Save the changes by clicking "Save changes"
- To test the function, click the "Test" tab and follow the instructions 
- Once the function is working, go to the "Publish" tab and click "Publish function". **Note:** if there already is 
  a function published that is live, it's possible to inspect the current live 
- To use the function a distribution has to be associated with it
- Click "Add association":
  - Select the relevant distribution (the associated domains should be listed with each distribution)
  - For "Event type", select "Viewer Request"
  - For "Cache behavior", select the behavior set up earlier that is associated with requests going to the S3 bucket, 
    in this case `discovery_chat_app/*`
  - Click "Add association"
- Click "View Distribution" to go back to the distribution page and wait for the distribution to be deployed.
- Once the distribution is deployed, try it out: [https://discovery.indra.bio/discovery_chat_app/](https://discovery.indra.bio/discovery_chat_app/)

#### Restrict access to bucket to only the CloudFront distribution (optional)

This can be done in order to restrict access to the bucket itself but allow cloudfront to serve the content.
